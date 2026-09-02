#!/usr/bin/env python3
"""
msh2Dto3D.py -- turn a 2-D gmsh mesh into a one-cell-thick 3-D mesh that
gmshToFoam and splitMeshRegions can actually use, and fix the OpenFOAM
boundary file afterwards.

WHY THIS EXISTS
---------------
gmshToFoam fed a 2-D .msh reports

    Cells: total:0  hex:0  prism:0  pyr:0  tet:0

and then dies. OpenFOAM needs VOLUME cells even for a planar case, so the mesh
has to be extruded by one layer with front/back marked `empty`. Doing that by
hand has cost this project a documented twelve distinct failures before the
first time step (SoEEDF/docs/setting-up-a-new-case.md), of which five belong to
this step alone. This script removes all five.

WHY IT EXTRUDES THE MESH AND NOT THE GEOMETRY
---------------------------------------------
Re-authoring the .geo to extrude the CAD is the obvious route and it is the one
that generates the traps:

  * under SetFactory("OpenCASCADE"), extruding two surfaces that SHARE a curve
    produces two COINCIDENT BUT SEPARATE lateral faces. The volumes are then
    non-conformal, splitMeshRegions cannot build the region interface, and both
    regions come back with `defaultFaces` and no coupling -- measured;
  * `Coherence;` renumbers entities and invalidates the Extrude return indices
    the physical groups were built from -- measured, face counts collapsed from
    103/160/338/352 to 3/1/7/184.

Both are CAD problems. Extruding the MESH has no CAD step, so neither can
happen: the 2-D mesh already shares nodes along the interface, and duplicating
nodes in z preserves that by construction.

ZERO IMPORTS BEYOND THE STANDARD LIBRARY
----------------------------------------
Deliberate. `python3` and `python` on PATH here resolve to an environment that
has neither `gmsh` nor `meshio`, so a script depending on either would fail for
the next person while working for whoever wrote it. MSH 2.2 ASCII is a trivial
format -- $PhysicalNames, $Nodes, $Elements, no geometry section -- so there is
nothing to gain by taking a dependency.

WHAT IT HANDLES FOR YOU
-----------------------
1. one layer, tri -> PRISM and quad -> hex. (The cells are prisms when the base
   mesh is triangular. gmsh's `Recombine` makes the LATERAL faces quads; it does
   not turn triangles into quads, and that surprises people.)
2. front/back physical groups, so the boundary file can mark them `empty`.
3. REGION INTERFACES ARE DETECTED AND EXCLUDED AUTOMATICALLY. An edge whose two
   adjacent faces lie in DIFFERENT physical surfaces is internal to the assembly
   and must NOT become a boundary patch -- splitMeshRegions builds
   gas_to_dielectric / dielectric_to_gas from INTERNAL faces between cellZones,
   so naming those faces makes the coupling impossible. This is inferred from
   the mesh; you do not have to know which curve tag it is. In the shipped
   needle mesh the group `air_dielectric` mixes 160 interface edges with 178
   genuine wall edges, so excluding the whole group would be wrong too.
4. physical groups are promoted a dimension: 2-D Physical Curve -> 3-D Physical
   Surface, 2-D Physical Surface -> 3-D Physical Volume. NAMES are preserved;
   TAGS renumber, because front/back are inserted. Everything downstream must
   therefore key off names, which changeDictionary does.
5. triangle winding is normalised so every prism is consistently oriented.

USAGE
-----
Convert:
    tools/msh2Dto3D.py mesh2D.msh -o mesh3D.msh --thickness 1e-5

Inspect before converting (recommended -- see --report below):
    tools/msh2Dto3D.py mesh2D.msh --report

Fix the OpenFOAM boundary file, AFTER gmshToFoam:
    tools/msh2Dto3D.py --fix-boundary constant/polyMesh/boundary

Full sequence for a two-region case:
    tools/msh2Dto3D.py mesh2D.msh -o mesh3D.msh
    gmshToFoam mesh3D.msh
    tools/msh2Dto3D.py --fix-boundary constant/polyMesh/boundary
    splitMeshRegions -cellZones -overwrite

ALWAYS RUN --report FIRST
-------------------------
An empty physical group is NOT an error in gmsh: it declares the group, exits 0,
and the mesh simply lacks those faces. That produced a mesh with 19 lateral
faces instead of thousands, and the only symptom was patches missing after
conversion. --report prints the element count per group so a zero is visible
before you spend a run finding it.
"""

import argparse
import sys
from collections import defaultdict

# MSH element type -> node count, for the types a 2-D mesh and its extrusion use.
NODES_PER_TYPE = {
    1: 2,    # 2-node line
    2: 3,    # 3-node triangle
    3: 4,    # 4-node quadrangle
    4: 4,    # 4-node tetrahedron
    5: 8,    # 8-node hexahedron
    6: 6,    # 6-node prism
    15: 1,   # 1-node point
}

LINE, TRI, QUAD, HEX, PRISM = 1, 2, 3, 5, 6


class Mesh:
    def __init__(self):
        self.names = {}          # (dim, tag) -> name
        self.nodes = {}          # id -> (x, y, z)
        self.elements = []       # (id, type, [tags], [node ids])


def read_msh(path):
    """Parse MSH 2.2 ASCII. Raises on anything it does not understand rather
    than guessing, because a silently mis-parsed mesh is worse than a stop."""
    m = Mesh()
    with open(path, "r") as fh:
        lines = fh.read().splitlines()

    i = 0
    n = len(lines)
    seen_format = False

    while i < n:
        s = lines[i].strip()

        if s == "$MeshFormat":
            ver = lines[i + 1].split()
            if not ver[0].startswith("2."):
                raise SystemExit(
                    "msh2Dto3D: this reader handles MSH 2.x ASCII; the file says"
                    " version %s.\n"
                    "  Re-export from gmsh with:  Mesh.MshFileVersion = 2.2;\n"
                    "  or on the CLI:             gmsh -3 -format msh22 ..."
                    % ver[0]
                )
            if len(ver) > 1 and ver[1] != "0":
                raise SystemExit(
                    "msh2Dto3D: the file is BINARY MSH (file-type %s). Re-export"
                    " as ASCII." % ver[1]
                )
            seen_format = True
            i += 2

        elif s == "$PhysicalNames":
            count = int(lines[i + 1])
            for k in range(count):
                parts = lines[i + 2 + k].split(None, 2)
                dim, tag = int(parts[0]), int(parts[1])
                m.names[(dim, tag)] = parts[2].strip().strip('"')
            i += 2 + count

        elif s == "$Nodes":
            count = int(lines[i + 1])
            for k in range(count):
                p = lines[i + 2 + k].split()
                m.nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
            i += 2 + count

        elif s == "$Elements":
            count = int(lines[i + 1])
            for k in range(count):
                p = lines[i + 2 + k].split()
                eid, etype, ntags = int(p[0]), int(p[1]), int(p[2])
                tags = [int(x) for x in p[3:3 + ntags]]
                nds = [int(x) for x in p[3 + ntags:]]
                if etype not in NODES_PER_TYPE:
                    raise SystemExit(
                        "msh2Dto3D: element %d has unsupported type %d."
                        % (eid, etype)
                    )
                if len(nds) != NODES_PER_TYPE[etype]:
                    raise SystemExit(
                        "msh2Dto3D: element %d (type %d) has %d nodes, expected"
                        " %d." % (eid, etype, len(nds), NODES_PER_TYPE[etype])
                    )
                m.elements.append((eid, etype, tags, nds))
            i += 2 + count

        else:
            i += 1

    if not seen_format:
        raise SystemExit("msh2Dto3D: no $MeshFormat section -- is this a .msh file?")
    if not m.nodes:
        raise SystemExit("msh2Dto3D: no nodes found.")
    if not m.elements:
        raise SystemExit("msh2Dto3D: no elements found.")
    return m


def phys_tag(tags):
    """Physical tag is the first tag in MSH 2.x. 0 means 'no physical group'."""
    return tags[0] if tags else 0


def classify_edges(mesh):
    """Split the 1-D elements into genuine boundary edges and region interfaces.

    An edge is an INTERFACE when its two adjacent 2-D faces belong to different
    physical surfaces. Such faces must stay INTERNAL so splitMeshRegions can
    build the region coupling from them.

    Returns (boundary_edges, interface_edges, edge_face_regions).
    """
    # For every undirected edge of every 2-D face, record the face's physical tag.
    edge_regions = defaultdict(set)
    for _, etype, tags, nds in mesh.elements:
        if etype not in (TRI, QUAD):
            continue
        p = phys_tag(tags)
        k = len(nds)
        for a in range(k):
            e = (nds[a], nds[(a + 1) % k])
            edge_regions[(min(e), max(e))].add(p)

    boundary, interface = [], []
    for el in mesh.elements:
        _, etype, _, nds = el
        if etype != LINE:
            continue
        key = (min(nds[0], nds[1]), max(nds[0], nds[1]))
        regions = edge_regions.get(key, set())
        if len(regions) > 1:
            interface.append(el)
        else:
            boundary.append(el)
    return boundary, interface, edge_regions


def signed_area(mesh, nds):
    """2-D signed area (xy) -- used only to normalise winding."""
    a = 0.0
    k = len(nds)
    for i in range(k):
        x1, y1, _ = mesh.nodes[nds[i]]
        x2, y2, _ = mesh.nodes[nds[(i + 1) % k]]
        a += x1 * y2 - x2 * y1
    return 0.5 * a


def report(mesh):
    boundary, interface, _ = classify_edges(mesh)

    per_group = defaultdict(int)
    for _, etype, tags, _ in mesh.elements:
        per_group[(1 if etype == LINE else 2, phys_tag(tags))] += 1

    iface_per_group = defaultdict(int)
    for _, _, tags, _ in interface:
        iface_per_group[phys_tag(tags)] += 1

    print("mesh: %d nodes, %d elements" % (len(mesh.nodes), len(mesh.elements)))
    print("")
    print("  %-22s %4s %8s %10s %s" % ("physical group", "dim", "elements",
                                       "interface", "becomes"))
    print("  " + "-" * 68)

    empty = []
    for (dim, tag), name in sorted(mesh.names.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        count = per_group.get((dim, tag), 0)
        nif = iface_per_group.get(tag, 0) if dim == 1 else 0
        if dim == 1:
            surviving = count - nif
            becomes = ("patch (%d faces)" % surviving) if surviving else "NOTHING"
        else:
            becomes = "cellZone -> region"
        print("  %-22s %4d %8d %10s %s"
              % (name, dim, count, (str(nif) if nif else "-"), becomes))
        if count == 0:
            empty.append(name)

    print("")
    print("  interface edges detected (excluded from patches): %d" % len(interface))
    print("  boundary edges -> lateral faces:                  %d" % len(boundary))

    if empty:
        print("")
        print("  *** WARNING: %d physical group(s) contain NO elements: %s"
              % (len(empty), ", ".join(empty)))
        print("      An empty group is not an error in gmsh -- it declares the")
        print("      group and exits 0 -- so this is normally a selection that")
        print("      matched nothing (e.g. `Surface In BoundingBox`). Those")
        print("      patches will be missing after gmshToFoam.")

    if not interface and sum(1 for (d, _) in mesh.names if d == 2) > 1:
        print("")
        print("  *** NOTE: more than one 2-D physical group, but NO interface")
        print("      edges were found. If those groups are meant to touch, the")
        print("      meshes are not conformal and splitMeshRegions will not be")
        print("      able to couple them.")


def convert(mesh, thickness, front_name, back_name):
    boundary, interface, _ = classify_edges(mesh)

    nmax = max(mesh.nodes)
    out_nodes = []
    for nid in sorted(mesh.nodes):
        x, y, z = mesh.nodes[nid]
        out_nodes.append((nid, x, y, z))
    for nid in sorted(mesh.nodes):
        x, y, z = mesh.nodes[nid]
        out_nodes.append((nid + nmax, x, y, z + thickness))

    def top(n):
        return n + nmax

    # --- physical groups: surfaces first, then volumes (matching gmsh's own
    # --- ordering habit), names preserved, tags renumbered.
    new_names = []          # (dim, newtag, name)
    tag_of = {}             # ('curve'|'surface', oldtag) -> newtag
    next_tag = 1

    surviving_curve_tags = sorted({phys_tag(t) for _, _, t, _ in boundary
                                   if phys_tag(t) != 0})
    for old in surviving_curve_tags:
        name = mesh.names.get((1, old), "curve_%d" % old)
        new_names.append((2, next_tag, name))
        tag_of[("curve", old)] = next_tag
        next_tag += 1

    front_tag = next_tag
    new_names.append((2, front_tag, front_name))
    next_tag += 1
    back_tag = next_tag
    new_names.append((2, back_tag, back_name))
    next_tag += 1

    surface_tags = sorted({phys_tag(t) for _, et, t, _ in mesh.elements
                           if et in (TRI, QUAD) and phys_tag(t) != 0})
    for old in surface_tags:
        name = mesh.names.get((2, old), "surface_%d" % old)
        new_names.append((3, next_tag, name))
        tag_of[("surface", old)] = next_tag
        next_tag += 1

    # --- elements
    out_elems = []
    eid = 1

    # cells: tri -> prism, quad -> hex, winding normalised so +z extrusion gives
    # a consistently oriented cell.
    for _, etype, tags, nds in mesh.elements:
        if etype not in (TRI, QUAD):
            continue
        p = tag_of.get(("surface", phys_tag(tags)))
        if p is None:
            continue                      # face in no physical group -> no cell
        v = list(nds)
        if signed_area(mesh, v) < 0:
            v = list(reversed(v))
        cell = v + [top(k) for k in v]
        out_elems.append((eid, PRISM if etype == TRI else HEX, p, p, cell))
        eid += 1

    # lateral faces from genuine boundary edges only
    for _, _, tags, nds in boundary:
        p = tag_of.get(("curve", phys_tag(tags)))
        if p is None:
            continue                      # edge in no physical group -> skip
        a, b = nds[0], nds[1]
        out_elems.append((eid, QUAD, p, p, [a, b, top(b), top(a)]))
        eid += 1

    # front (z = 0) and back (z = thickness)
    for _, etype, tags, nds in mesh.elements:
        if etype not in (TRI, QUAD):
            continue
        if tag_of.get(("surface", phys_tag(tags))) is None:
            continue
        v = list(nds)
        if signed_area(mesh, v) < 0:
            v = list(reversed(v))
        out_elems.append((eid, etype, front_tag, front_tag, v))
        eid += 1
        out_elems.append((eid, etype, back_tag, back_tag,
                          [top(k) for k in v]))
        eid += 1

    return out_nodes, out_elems, new_names, len(interface)


def write_msh(path, out_nodes, out_elems, new_names):
    with open(path, "w") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write("$PhysicalNames\n%d\n" % len(new_names))
        for dim, tag, name in new_names:
            fh.write('%d %d "%s"\n' % (dim, tag, name))
        fh.write("$EndPhysicalNames\n")
        fh.write("$Nodes\n%d\n" % len(out_nodes))
        for nid, x, y, z in out_nodes:
            fh.write("%d %.17g %.17g %.17g\n" % (nid, x, y, z))
        fh.write("$EndNodes\n")
        fh.write("$Elements\n%d\n" % len(out_elems))
        for eid, etype, phys, elem, nds in out_elems:
            fh.write("%d %d 2 %d %d %s\n"
                     % (eid, etype, phys, elem,
                        " ".join(str(k) for k in nds)))
        fh.write("$EndElements\n")


def fix_boundary(path, empty_names):
    """Post-gmshToFoam fixes to constant/polyMesh/boundary.

    gmsh cannot express either of these, so they have to happen on the OpenFOAM
    side:

      * `physicalType patch;` -- gmshToFoam writes it and changeDictionary
        chokes on it;
      * front/back must be `empty`, which is what makes this a standard
        OpenFOAM 2-D planar case. gmshToFoam leaves them as `patch`.
    """
    with open(path, "r") as fh:
        lines = fh.read().splitlines()

    out = []
    dropped = 0
    retyped = []
    pending = None      # patch name whose block we are inside

    for idx, line in enumerate(lines):
        stripped = line.strip()

        if stripped.startswith("physicalType"):
            dropped += 1
            continue

        # A bare word on its own line, followed by '{', is a patch name.
        if (stripped and " " not in stripped and not stripped.endswith(";")
                and stripped not in ("{", "}", "(", ")")
                and idx + 1 < len(lines) and lines[idx + 1].strip() == "{"):
            pending = stripped

        if pending in empty_names and stripped.startswith("type"):
            indent = line[:len(line) - len(line.lstrip())]
            out.append("%stype            empty;" % indent)
            retyped.append(pending)
            pending = None
            continue

        if stripped == "}":
            pending = None

        out.append(line)

    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")

    print("fix-boundary: %s" % path)
    print("  removed %d `physicalType` entr%s"
          % (dropped, "y" if dropped == 1 else "ies"))
    if retyped:
        print("  set `type empty;` on: %s" % ", ".join(retyped))
    for name in empty_names:
        if name not in retyped:
            print("  NOTE: no patch named `%s` in this file -- nothing to do"
                  % name)


def main():
    ap = argparse.ArgumentParser(
        description="Extrude a 2-D gmsh mesh into a one-cell-thick 3-D mesh for "
                    "OpenFOAM, and fix the boundary file afterwards.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run --report first: an empty physical group is not an error in "
               "gmsh, and the only symptom is patches missing after gmshToFoam.")
    ap.add_argument("input", nargs="?", help="the 2-D .msh file (MSH 2.x ASCII)")
    ap.add_argument("-o", "--output", help="the 3-D .msh to write")
    ap.add_argument("-t", "--thickness", type=float, default=1e-5,
                    help="extrusion thickness [m]; arbitrary for `empty` "
                         "patches (default 1e-5)")
    ap.add_argument("--front-name", default="front",
                    help="physical name for the z=0 face (default front)")
    ap.add_argument("--back-name", default="back",
                    help="physical name for the z=thickness face (default back)")
    ap.add_argument("--report", action="store_true",
                    help="print per-group element counts and detected "
                         "interfaces, then exit without writing")
    ap.add_argument("--fix-boundary", metavar="BOUNDARY",
                    help="post-gmshToFoam: strip `physicalType` and set "
                         "front/back to `empty` in constant/polyMesh/boundary")
    args = ap.parse_args()

    if args.fix_boundary:
        fix_boundary(args.fix_boundary, {args.front_name, args.back_name})
        return 0

    if not args.input:
        ap.error("an input .msh is required (or use --fix-boundary)")

    mesh = read_msh(args.input)

    if args.report:
        report(mesh)
        return 0

    if not args.output:
        ap.error("-o/--output is required when converting")

    if args.thickness <= 0:
        ap.error("--thickness must be positive")

    out_nodes, out_elems, new_names, n_iface = convert(
        mesh, args.thickness, args.front_name, args.back_name)
    write_msh(args.output, out_nodes, out_elems, new_names)

    ncells = sum(1 for e in out_elems if e[1] in (PRISM, HEX))
    print("wrote %s" % args.output)
    print("  %d nodes, %d elements (%d cells)"
          % (len(out_nodes), len(out_elems), ncells))
    print("  thickness %g m, front=`%s` back=`%s`"
          % (args.thickness, args.front_name, args.back_name))
    print("  %d interface edges kept INTERNAL (not made into patches)" % n_iface)
    print("")
    print("  physical groups written:")
    for dim, tag, name in new_names:
        kind = "volume" if dim == 3 else "surface"
        print("    %d %-8s %s" % (tag, kind, name))
    print("")
    print("  next:")
    print("    gmshToFoam %s" % args.output)
    print("    %s --fix-boundary constant/polyMesh/boundary" % sys.argv[0])
    print("    splitMeshRegions -cellZones -overwrite")
    return 0


if __name__ == "__main__":
    sys.exit(main())
