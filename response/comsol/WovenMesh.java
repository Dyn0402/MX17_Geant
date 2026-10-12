/*
 * WovenMesh.java -- electrostatic field of one cell of the TRUE plain-weave micromesh.
 *
 * The production map (meshcell/solve_fieldmap.py) stands the weave in by two flat
 * layers of straight crossed wires at z = +-9 um. Here every wire is a sinusoid,
 * z = +-A sin(pi s / p), swept with its real 9.5 um radius, over and under its
 * crossing partners, with the same 1 um kiss overlap (A = r - overlap/2 = 9.0 um).
 *
 * Cell: x, y in [-p/2, p/2], wires along the four side faces (half-wires). The planes
 * through wire centres are exact mirror planes of a plain weave, so this cell with
 * zero-charge (Neumann) side walls IS the infinite weave, and Garfield must load the
 * map with mirror periodicity (EnableMirrorPeriodicityX/Y), not plain periodicity:
 * the weave's translation period is 2p.
 *
 *   x-wire at y = -p/2:  z = -A sin(pi x/p)      x-wire at y = +p/2:  z = +A sin(pi x/p)
 *   y-wire at x = -p/2:  z = +A sin(pi y/p)      y-wire at x = +p/2:  z = -A sin(pi y/p)
 *
 * so at every crossing the two wires sit at +-A, as in a weave.
 *
 * Boundaries as in solve_fieldmap.py: mesh 0 V, ESL anode z = -168.5 um at +V_MESH
 * (150 um pillars to the mesh underside, -18.5 um), drift slice top at +218.5 um. The
 * top carries a surface charge that makes the far drift field exactly E_DRIFT (this
 * replaces solve_fieldmap.py's two-solve superposition and is exact for a linear problem).
 *
 * Environment (all optional):
 *   V_MESH    anode potential [V]                   default 490
 *   E_DRIFT   drift field [V/cm]                    default 333
 *   LC_WIRE   max element size on wires [um]        default 1.0
 *   LC_FAR    max element size far away [um]        default 8
 *   MAP       1 = print the Garfield grid as "CSV," lines (run_lxplus.sh writes them to
 *             $OUT; COMSOL's sandbox forbids file I/O from model code, both ways), then
 *             woven_map.py assemble turns them into the ComponentGrid file
 */
import com.comsol.model.*;
import com.comsol.model.util.*;
import com.comsol.model.physics.*;
import java.util.*;

public class WovenMesh {
  static String env(String k, String d) { String v = System.getenv(k); return (v == null || v.isEmpty()) ? d : v; }

  public static Model run() throws Exception {
    double vMesh = Double.parseDouble(env("V_MESH", "490"));
    double eDrift = Double.parseDouble(env("E_DRIFT", "333"));
    double lcWire = Double.parseDouble(env("LC_WIRE", "1.0"));
    double lcFar = Double.parseDouble(env("LC_FAR", "8"));

    final double P = 67.0, RW = 9.5, OVL = 1.0, A = RW - OVL / 2;   // A = 9.0
    final double ZU = -(A + RW), ZA = ZU - 150.0, ZT = (A + RW) + 200.0;   // -18.5, -168.5, 218.5

    Model m = ModelUtil.create("Model");
    m.param().set("p", P + "[um]");
    m.param().set("A", A + "[um]");
    ModelNode comp = m.component().create("comp1", true);
    GeomSequence g = comp.geom().create("geom1", 3);
    g.lengthUnit("um");

    // sign of each wire's sinusoid (see header) and its fixed transverse coordinate
    double[] xwY = {-P / 2, P / 2}, xwS = {-1, 1};    // x-wires: y = xwY, z = xwS*A sin(pi x/p)
    double[] ywX = {-P / 2, P / 2}, ywS = {1, -1};    // y-wires: x = ywX, z = ywS*A sin(pi y/p)
    double s0 = -1.5 * P, s1 = 1.5 * P;               // sin' = 0 at both ends: the section is planar
    String[] wires = new String[4];
    for (int j = 0; j < 2; j++) {
      String t = "xw" + j;
      double z0 = xwS[j] * A * Math.sin(Math.PI * s0 / P);
      GeomFeature wp = g.create("wp_" + t, "WorkPlane");
      wp.set("quickplane", "yz");
      wp.set("quickx", s0);
      wp.geom().create("c", "Circle").set("pos", new double[]{xwY[j], z0});
      wp.geom().feature("c").set("r", RW);
      GeomFeature pc = g.create("pc_" + t, "ParametricCurve");
      pc.set("parmin", s0); pc.set("parmax", s1);
      pc.set("coord", new String[]{"s", String.valueOf(xwY[j]), xwS[j] + "*" + A + "*sin(pi*s/" + P + ")"});
      wires[j] = sweep(g, "sw_" + t, "wp_" + t, "pc_" + t);
    }
    for (int i = 0; i < 2; i++) {
      String t = "yw" + i;
      double z0 = ywS[i] * A * Math.sin(Math.PI * s0 / P);
      GeomFeature wp = g.create("wp_" + t, "WorkPlane");
      wp.set("quickplane", "zx");
      wp.set("quicky", s0);
      wp.geom().create("c", "Circle").set("pos", new double[]{z0, ywX[i]});   // zx plane: (z, x)
      wp.geom().feature("c").set("r", RW);
      GeomFeature pc = g.create("pc_" + t, "ParametricCurve");
      pc.set("parmin", s0); pc.set("parmax", s1);
      pc.set("coord", new String[]{String.valueOf(ywX[i]), "s", ywS[i] + "*" + A + "*sin(pi*s/" + P + ")"});
      wires[2 + i] = sweep(g, "sw_" + t, "wp_" + t, "pc_" + t);
    }
    GeomFeature blk = g.create("blk", "Block");
    blk.set("pos", new double[]{-P / 2, -P / 2, ZA});
    blk.set("size", new double[]{P, P, ZT - ZA});
    GeomFeature dif = g.create("dif", "Difference");
    dif.selection("input").set("blk");
    dif.selection("input2").set(wires);
    g.run();
    System.out.println("[woven] geometry: " + g.getNDomains() + " domains, " + g.getNBoundaries() + " boundaries");

    SelectionFeature wsel = comp.selection().create("s_wire", "Box");
    wsel.set("entitydim", 2);
    wsel.set("xmin", -P); wsel.set("xmax", P); wsel.set("ymin", -P); wsel.set("ymax", P);
    wsel.set("zmin", ZU - 0.5); wsel.set("zmax", -ZU + 0.5);
    wsel.set("condition", "inside");
    SelectionFeature bot = comp.selection().create("s_bot", "Box");
    bot.set("entitydim", 2);
    bot.set("xmin", -P); bot.set("xmax", P); bot.set("ymin", -P); bot.set("ymax", P);
    bot.set("zmin", ZA - 0.01); bot.set("zmax", ZA + 0.01);
    bot.set("condition", "inside");
    SelectionFeature top = comp.selection().create("s_top", "Box");
    top.set("entitydim", 2);
    top.set("xmin", -P); top.set("xmax", P); top.set("ymin", -P); top.set("ymax", P);
    top.set("zmin", ZT - 0.01); top.set("zmax", ZT + 0.01);
    top.set("condition", "inside");

    Physics es = comp.physics().create("es", "Electrostatics", "geom1");
    es.create("gnd", "Ground", 2).selection().named("s_wire");
    PhysicsFeature an = es.create("anode", "ElectricPotential", 2);
    an.selection().named("s_bot");
    an.set("V0", vMesh + "[V]");
    // exterior boundary: -n.D = rho_s with n = +z, so D_z = eps0*E_drift needs rho_s = -eps0*E_drift
    PhysicsFeature sc = es.create("drift", "SurfaceChargeDensity", 2);
    sc.selection().named("s_top");
    sc.set("rhoqs", "-epsilon0_const*" + eDrift + "[V/cm]");

    MeshSequence ms = comp.mesh().create("mesh1");
    ms.create("sz", "Size");
    ms.feature("sz").set("custom", "on");
    ms.feature("sz").set("hmax", lcFar + "[um]");
    ms.feature("sz").set("hmaxactive", true);
    ms.feature("sz").set("hmin", "0.1[um]");
    ms.feature("sz").set("hminactive", true);
    ms.feature("sz").set("hgrad", 1.3);
    ms.feature("sz").set("hgradactive", true);
    ms.create("ftet", "FreeTet");
    ms.feature("ftet").create("sw", "Size");
    ms.feature("ftet").feature("sw").selection().geom("geom1", 2);
    ms.feature("ftet").feature("sw").selection().named("s_wire");
    ms.feature("ftet").feature("sw").set("custom", "on");
    ms.feature("ftet").feature("sw").set("hmax", lcWire + "[um]");
    ms.feature("ftet").feature("sw").set("hmaxactive", true);
    ms.run();
    System.out.println("[woven] mesh elements: " + comp.mesh("mesh1").getNumElem());

    m.study().create("std1");
    m.study("std1").create("stat", "Stationary");
    m.study("std1").run();

    // gates (cf. solve_fieldmap.py S1/S2/S3), printed for the log
    gate(m, "S1 Ez drift 40um below top [V/cm]", "es.Ez", 0.25 * P, 0.25 * P, ZT - 40);
    gate(m, "S1 Ez drift 40um below top, centre [V/cm]", "es.Ez", 0, 0, ZT - 40);
    gate(m, "S2 Ez amp bulk z=-100 [V/cm]", "es.Ez", 0, 0, -100);
    gate(m, "S2 Ez amp bulk z=-100 corner [V/cm]", "es.Ez", 0.4 * P, 0.4 * P, -100);
    gate(m, "V at drift top [V]", "V", 0, 0, ZT);
    // the potential just off a wire: surface point of the x-wire at y=-p/2, x=0 (z=0), 0.4 um out
    gate(m, "S3 V 0.4um off x-wire [V]", "V", 0, -P / 2 + RW + 0.4, 0);

    if (env("MAP", "0").equals("1")) writeMap(m, P, RW, A);
    return m;
  }

  /*
   * The production grid of meshcell/solve_fieldmap.py: 1 um step, x, y in [-p/2, p/2],
   * z in [-168, 218], z outer, x, y inner. Flags as there (_sample_plane): 1 gas,
   * 2 in-wire within one step of the surface (evaluated 0.4 um outside the nearest
   * wire surface, along the radius), 0 deeper. woven_map.py classify() is the same
   * logic in numpy and cross-checks the counts.
   */
  static void writeMap(Model m, double P, double RW, double A) {
    final double STEP = 1.0, STANDOFF = 0.4;
    int nx = (int) Math.round(P / STEP) + 1, nz = 387;
    double[][] wires = {{0, -P / 2, -1}, {0, P / 2, 1}, {1, -P / 2, 1}, {1, P / 2, -1}};  // axis(0=x,1=y), c, sign
    NumericalFeature ip = m.result().numerical().create("mapint", "Interp");
    ip.set("data", "dset1");
    ip.set("expr", new String[]{"V", "es.Ex", "es.Ey", "es.Ez"});
    ip.set("unit", new String[]{"V", "V/cm", "V/cm", "V/cm"});
    System.out.println("CSV,# WovenMesh grid: ix iy iz flag V Ex Ey Ez (V, V/cm); x = -p/2 + ix um, z = -168 + iz um");
    int n = nx * nx;
    for (int iz = 0; iz < nz; iz++) {
      double z = -168.0 + iz * STEP;
      int[] flag = new int[n];
      double[][] q = new double[3][n];
      int ne = 0;
      int[] idx = new int[n];
      for (int ix = 0; ix < nx; ix++) for (int iy = 0; iy < nx; iy++) {
        double x = -P / 2 + ix * STEP, y = -P / 2 + iy * STEP;
        double dmin = 1e9; double[] cbest = null;
        for (double[] w : wires) {
          boolean alongX = w[0] == 0;
          double s0 = alongX ? x : y, t = alongX ? y : x, k = Math.PI / P, sg = w[2], s = s0;
          for (int it = 0; it < 8; it++) {
            double zs = sg * A * Math.sin(k * s), dz = sg * A * k * Math.cos(k * s), d2z = -sg * A * k * k * Math.sin(k * s);
            s -= ((s - s0) + (zs - z) * dz) / (1 + dz * dz + (zs - z) * d2z);
          }
          double zs = sg * A * Math.sin(k * s);
          double d = Math.sqrt((s - s0) * (s - s0) + (t - w[1]) * (t - w[1]) + (z - zs) * (z - zs));
          if (d < dmin) { dmin = d; cbest = alongX ? new double[]{s, w[1], zs} : new double[]{w[1], s, zs}; }
        }
        int g = ix * nx + iy;
        flag[g] = dmin >= RW ? 1 : (dmin >= RW - STEP ? 2 : 0);
        if (flag[g] == 0) continue;
        double sc = flag[g] == 2 ? (RW + STANDOFF) / Math.max(dmin, 1e-6) : 1.0;
        q[0][ne] = flag[g] == 2 ? cbest[0] + (x - cbest[0]) * sc : x;
        q[1][ne] = flag[g] == 2 ? cbest[1] + (y - cbest[1]) * sc : y;
        q[2][ne] = flag[g] == 2 ? cbest[2] + (z - cbest[2]) * sc : z;
        idx[ne++] = g;
      }
      double[][] c = {Arrays.copyOf(q[0], ne), Arrays.copyOf(q[1], ne), Arrays.copyOf(q[2], ne)};
      ip.set("coord", c);
      double[][] v = ip.getReal();
      if (iz == 0) System.out.println("[woven] Interp getReal shape " + v.length + " x " + v[0].length + " for " + ne + " points");
      // Interp with a coord matrix returns [point][expr]
      if (v.length != ne || v[0].length != 4) throw new RuntimeException("Interp returned " + v.length + " x " + v[0].length + ", expected " + ne + " x 4");
      StringBuilder sb = new StringBuilder();
      int e = 0;
      for (int g = 0; g < n; g++) {
        sb.append("CSV,").append(g / nx).append(' ').append(g % nx).append(' ').append(iz).append(' ').append(flag[g]);
        if (flag[g] > 0) { for (int r = 0; r < 4; r++) sb.append(' ').append((float) v[e][r]); e++; }
        else sb.append(" 0 0 0 0");
        sb.append('\n');
      }
      System.out.print(sb);
    }
    System.out.println("[woven] map written");
  }

  static String sweep(GeomSequence g, String tag, String wp, String pc) {
    g.run(wp); g.run(pc);       // a selection can only name objects that already exist
    GeomFeature sw = g.create(tag, "Sweep");
    sw.selection("face").set(wp, 1);
    sw.selection("edge").set(pc, 1);
    sw.set("keep", false);
    return tag;
  }

  static void gate(Model m, String label, String expr, double x, double y, double z) {
    String tag = "pt" + Math.abs(label.hashCode());
    NumericalFeature ip = m.result().numerical().create(tag, "Interp");
    ip.set("data", "dset1");
    ip.set("expr", new String[]{expr});
    ip.set("unit", new String[]{expr.equals("V") ? "V" : "V/cm"});
    ip.set("coord", new double[][]{{x}, {y}, {z}});
    double[][] v = ip.getReal();
    System.out.println("[woven] GATE " + label + " = " + v[0][0]);
  }

  public static void main(String[] args) throws Exception { run(); }
}
