/*
 * Induce3D.java -- direct time-domain induction in 3D: a point avalanche on the real
 * resistive STRIPS over the real checkerboard pad comb.
 *
 * Induce2D had the strips but as infinite floating lines (no spreading); InduceAxi had
 * spreading but over a uniform sheet with ring pads. This one has both: the landed charge
 * diffuses ALONG its 550 um strip (strips on an 800 um pitch run along y, isolated from
 * their neighbours by 250 um gaps), and it is read by the comb channels the board has
 * (response/common/channel_map.py, kernels.py):
 *
 *   a Y channel = one pad ROW j,    pads at columns i with (i + j) EVEN
 *   an X channel = one pad COLUMN i, pads at rows    j with (i + j) ODD
 *
 * Pads 680 um on a 780 um pitch, 100 um inter-pad trenches filled with glue (floating, as
 * W2/V6 require). Same stack, method and source as Induce2D/InduceAxi: the avalanche is a
 * prescribed convection current, COMSOL solves div(sigma E + eps dE/dt + Je) = 0, and the
 * channel Terminals at 0 V carry the induced currents. No weighting-field theorem.
 *
 * Quarter model. The avalanche sits on the axis x = y = 0, over the centre of strip 0 and
 * of pad (0,0) (a Y pad). x = 0 and y = 0 are exact mirror planes of the strips, of the
 * pads and of the checkerboard (mirroring keeps the parity of i + j), so the model holds
 * x, y >= 0 with zero-flux symmetry walls, and the full-detector channel currents are
 *   X0 = 4 x quarter, X_i>0 = 2 x quarter (its y mirror; x mirror is channel X_-i),
 *   Y0 = 4 x quarter, Y_j>0 = 2 x quarter,
 * applied before printing. The source term uses the full Q, of which the quarter holds Q/4.
 *
 * Outer walls x = XW (mid-gap after pad column 3) and y = YW (mid-gap after pad row 5,
 * 4.3 mm along the strip) are zero-flux: within the solve window the charge spreads
 * ~1.5 mm along the strip, so the strip ends (bus bars, far away on the real board) are
 * never reached.
 *
 * Geometry is drawn once as a 2D partition (strip edges, pad edges, a refinement square at
 * the avalanche) and extruded through every layer, so every domain is a prism and the mesh
 * is a triangle mesh swept through the stack.
 *
 * Environment (all optional):
 *   RHO_S     sheet resistivity [ohm/sq]          default 2e6
 *   CASE      pair | deposit                       default pair
 *   T_ION_NS  ion transit across the full gap      default 340
 *   T_END_NS  end of the solve                     default 1500
 *   VE_UM_NS  electron velocity [um/ns]            default 15 (use 150 for deposit)
 *   ORDER     shape order of the potential         default 2
 *   H_FAR_UM, H_NEAR_UM, NZ_GAS                    lateral mesh far/near, gas layers (120, 8, 14)
 *
 * Output: CSV lines on stdout prefixed "CSV,"; run_lxplus.sh writes them to $OUT. Columns are
 * full-detector channel currents X0..X3, Y0..Y5 in A for Q = 1e-15 C.
 */
import com.comsol.model.*;
import com.comsol.model.util.*;
import com.comsol.model.physics.*;
import java.util.*;

public class Induce3D {
  static String env(String k, String d) { String v = System.getenv(k); return (v == null || v.isEmpty()) ? d : v; }

  public static Model run() throws Exception {
    double rhoS = Double.parseDouble(env("RHO_S", "2e6"));
    String kase = env("CASE", "pair");
    double tIon = Double.parseDouble(env("T_ION_NS", "340"));
    double tEnd = Double.parseDouble(env("T_END_NS", "1500"));
    double ve = Double.parseDouble(env("VE_UM_NS", "15"));
    int order = Integer.parseInt(env("ORDER", "2"));
    double hFar = Double.parseDouble(env("H_FAR_UM", "120"));
    double hNear = Double.parseDouble(env("H_NEAR_UM", "8"));
    int nzGas = Integer.parseInt(env("NZ_GAS", "14"));

    // stack [um] -- response/common/constants.py, solver/v6_pad_gaps.py
    final double GAP = 150, T_ESL = 10, T_KAP = 50, T_GLUE = 18.76, T_CU = 26, T_FR4 = 264;
    final double PAD_P = 780, PAD_W = 680, ESL_P = 800, ESL_W = 550;
    final int NCOL = 4, NROW = 6;                          // pad columns i = 0..3, rows j = 0..5
    final double XW = (NCOL - 1) * PAD_P + PAD_P / 2;      // 2730: mid-gap after column 3
    final double YW = (NROW - 1) * PAD_P + PAD_P / 2;      // 4290: mid-gap after row 5
    final double RREF = 150;                               // refinement square at the avalanche
    // z levels, bottom up
    final double zB = -(T_ESL + T_KAP + T_GLUE + T_CU + T_FR4);
    final double zT = zB + T_FR4, zP = zT + T_CU, zG = zP + T_GLUE, zK = zG + T_KAP;   // zK + T_ESL = 0

    Model m = ModelUtil.create("Model");
    m.param().set("rho_s", rhoS + "[ohm]");
    m.param().set("t_esl", T_ESL + "[um]");
    m.param().set("gap", GAP + "[um]");
    m.param().set("Lz", "14.8[um]", "ion-birth decay length 1/alpha (ions.py)");
    m.param().set("Lb", "10[um]", "deposit case: electrons born within this of the mesh");
    m.param().set("sig0", "33[um]", "avalanche lateral sigma (T7)");
    m.param().set("vi", "gap/" + tIon + "[ns]");
    m.param().set("ve", ve + "[um/ns]", "electron velocity");
    m.param().set("Q", "1e-15[C]", "avalanche charge (full detector); output is per this Q");
    m.param().set("ion_on", kase.equals("deposit") ? "0" : "1");
    m.param().set("dep", kase.equals("deposit") ? "1" : "0");

    ModelNode comp = m.component().create("comp1", true);
    GeomSequence g = comp.geom().create("geom1", 3);
    g.lengthUnit("um");
    GeomFeature wp = g.create("wp1", "WorkPlane");
    wp.set("quickplane", "xy");
    wp.set("quickz", zB);
    GeomSequence w = wp.geom();
    rect2(w, "dom", 0, 0, XW, YW);
    rect2(w, "ref", 0, 0, RREF, RREF);
    List<double[]> strips = new ArrayList<>();
    for (int k = 0; k * ESL_P - ESL_W / 2 < XW; k++) {
      double a = Math.max(0, k * ESL_P - ESL_W / 2), b = Math.min(XW, k * ESL_P + ESL_W / 2);
      strips.add(new double[]{a, b});
      rect2(w, "st" + k, a, 0, b - a, YW);
    }
    for (int i = 0; i < NCOL; i++) for (int j = 0; j < NROW; j++) {
      double[] p = pad(i, j, PAD_P, PAD_W);
      rect2(w, "pd" + i + "_" + j, p[0], p[2], p[1] - p[0], p[3] - p[2]);
    }
    g.run("wp1");
    GeomFeature ext = g.create("ext1", "Extrude");
    ext.selection("input").set("wp1");
    ext.set("distance", new double[]{T_FR4, T_FR4 + T_CU, T_FR4 + T_CU + T_GLUE,
                                     T_FR4 + T_CU + T_GLUE + T_KAP, T_FR4 + T_CU + T_GLUE + T_KAP + T_ESL,
                                     T_FR4 + T_CU + T_GLUE + T_KAP + T_ESL + GAP});
    g.run();
    System.out.println("[induce3d] geometry: " + g.getNDomains() + " domains, " + g.getNBoundaries() + " boundaries");

    // selections
    double e = 0.01, big = 1e5;
    box(comp, "s_gas", 3, -e, big, -e, big, -e, GAP + e);
    box(comp, "s_esllay", 3, -e, big, -e, big, zK - e, e);
    box(comp, "s_kap", 3, -e, big, -e, big, zG - e, zK + e);
    box(comp, "s_glue", 3, -e, big, -e, big, zP - e, zG + e);
    box(comp, "s_culay", 3, -e, big, -e, big, zT - e, zP + e);
    box(comp, "s_fr4", 3, -e, big, -e, big, zB - e, zT + e);
    box(comp, "s_mesh", 2, -e, big, -e, big, GAP - e, GAP + e);
    box(comp, "s_bot", 2, -e, big, -e, big, zB - e, zB + e);
    List<String> stTags = new ArrayList<>();
    for (int k = 0; k < strips.size(); k++) {
      double[] s = strips.get(k);
      box(comp, "s_st" + k, 3, s[0] - e, s[1] + e, -e, big, zK - e, e);
      stTags.add("s_st" + k);
    }
    union(comp, "s_esl", 3, stTags);
    List<String> padTags = new ArrayList<>();
    Map<String, List<String>> chan = new LinkedHashMap<>();
    for (int i = 0; i < NCOL; i++) chan.put("X" + i, new ArrayList<>());
    for (int j = 0; j < NROW; j++) chan.put("Y" + j, new ArrayList<>());
    for (int i = 0; i < NCOL; i++) for (int j = 0; j < NROW; j++) {
      double[] p = pad(i, j, PAD_P, PAD_W);
      String t = "s_pd" + i + "_" + j;
      box(comp, t, 3, p[0] - e, p[1] + e, p[2] - e, p[3] + e, zT - e, zP + e);
      padTags.add(t);
      chan.get(((i + j) % 2 == 0) ? "Y" + j : "X" + i).add(t);
    }
    union(comp, "s_pads", 3, padTags);
    SelectionFeature all = comp.selection().create("s_all", "Explicit");
    all.all();
    diff(comp, "s_phys", 3, "s_all", "s_pads");
    diff(comp, "s_trench", 3, "s_culay", "s_pads");

    Physics ec = comp.physics().create("ec", "ConductiveMedia", "geom1");
    ec.selection().named("s_phys");
    ec.prop("ShapeProperty").set("order_electricpotential", order);
    ec.feature("cucns1").set("sigma_mat", "userdef");
    ec.feature("cucns1").set("sigma", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    ec.feature("cucns1").set("epsilonr_mat", "userdef");
    ec.feature("cucns1").set("epsilonr", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    medium(ec, "cu_esl", "s_esl", "1/(rho_s*t_esl)", "3");      // gaps between strips stay gas
    medium(ec, "cu_kap", "s_kap", "0", "3.5");
    medium(ec, "cu_glue", "s_glue", "0", "3.2");
    medium(ec, "cu_tr", "s_trench", "0", "3.2");
    medium(ec, "cu_fr4", "s_fr4", "0", "4.4");
    union(comp, "s_gnd", 2, Arrays.asList("s_mesh", "s_bot"));
    ec.create("gnd1", "Ground", 2).selection().named("s_gnd");
    List<String> names = new ArrayList<>(chan.keySet());
    for (String c : names) {
      SelectionFeature adj = comp.selection().create("s_b" + c, "Adjacent");
      adj.set("input", chan.get(c).toArray(new String[0]));
      PhysicsFeature t = ec.create("term_" + c, "Terminal", 2);
      t.selection().named("s_b" + c);
      t.set("TerminalType", "Voltage");
      t.set("V0", "0");
      t.set("TerminalName", c);
    }

    // the avalanche as a convection current, vertical only; abs() keeps exponents bounded
    m.func().create("pz", "Analytic");
    m.func("pz").set("funcname", "pz");
    m.func("pz").set("expr", "(s>=0)*(s<=gap)*exp(-abs(s)/Lz)/(Lz*(1-exp(-gap/Lz)))");
    m.func("pz").set("args", new String[]{"s"});
    m.func("pz").set("argunit", "m");
    m.func("pz").set("fununit", "1/m");
    m.func().create("pb", "Analytic");
    m.func("pb").set("funcname", "pb");
    m.func("pb").set("expr", "(s>=0)*(s<=gap)*exp(-abs(gap-s)/Lb)/(Lb*(1-exp(-gap/Lb)))");
    m.func("pb").set("args", new String[]{"s"});
    m.func("pb").set("argunit", "m");
    m.func("pb").set("fununit", "1/m");
    comp.variable().create("var1");
    comp.variable("var1").selection().geom("geom1", 3);
    comp.variable("var1").selection().named("s_gas");
    comp.variable("var1").set("gr", "exp(-(x^2+y^2)/(2*sig0^2))/(2*pi*sig0^2)");
    comp.variable("var1").set("Jz", "Q*gr*(ion_on*vi*pz(z-vi*t) + ve*((1-dep)*pz(z+ve*t) + dep*pb(z+ve*t)))");
    PhysicsFeature je = ec.create("ecd1", "ExternalCurrentDensity", 3);
    je.selection().named("s_gas");
    je.set("Je", new String[]{"0", "0", "Jz"});

    // mesh: triangles on the bottom face, refined at the avalanche, swept up the stack
    MeshSequence ms = comp.mesh().create("mesh1");
    box(comp, "s_botf", 2, -e, big, -e, big, zB - e, zB + e);
    box(comp, "s_reff", 2, -e, RREF + e, -e, RREF + e, zB - e, zB + e);
    ms.create("ftri", "FreeTri");
    ms.feature("ftri").selection().named("s_botf");
    ms.feature("ftri").create("sz", "Size");
    ms.feature("ftri").feature("sz").set("custom", "on");
    ms.feature("ftri").feature("sz").set("hmax", hFar + "[um]");
    ms.feature("ftri").feature("sz").set("hmaxactive", true);
    ms.feature("ftri").feature("sz").set("hgrad", 1.25);
    ms.feature("ftri").feature("sz").set("hgradactive", true);
    ms.feature("ftri").create("szr", "Size");
    ms.feature("ftri").feature("szr").selection().named("s_reff");
    ms.feature("ftri").feature("szr").set("custom", "on");
    ms.feature("ftri").feature("szr").set("hmax", hNear + "[um]");
    ms.feature("ftri").feature("szr").set("hmaxactive", true);
    ms.create("swe", "Sweep");
    sweepLayers(ms, "dg", "s_gas", nzGas);
    sweepLayers(ms, "de", "s_esllay", 2);
    sweepLayers(ms, "dk", "s_kap", 3);
    sweepLayers(ms, "dl", "s_glue", 1);
    sweepLayers(ms, "dc", "s_culay", 1);
    sweepLayers(ms, "df", "s_fr4", 4);
    ms.run();
    System.out.println("[induce3d] mesh elements: " + comp.mesh("mesh1").getNumElem());

    m.study().create("std1");
    m.study("std1").create("time", "Transient");
    m.study("std1").feature("time").set("tlist", "range(0,0.1[ns],2[ns]) range(3[ns],1[ns]," + tEnd + "[ns])");
    m.study("std1").createAutoSequences("all");
    SolverFeature td = m.sol("sol1").feature("t1");
    td.set("tstepsbdf", "free");
    td.set("maxstepconstraintbdf", "const");
    td.set("maxstepbdf", "1[ns]");
    td.set("initialstepbdfactive", true);
    td.set("initialstepbdf", "1e-13");
    td.set("rtol", "1e-4");
    m.sol("sol1").feature("v1").set("scalemethod", "manual");
    m.sol("sol1").feature("v1").set("scaleval", "1e-1");
    m.sol("sol1").runAll();

    NumericalFeature ev = m.result().numerical().create("gev1", "EvalGlobal");
    String[] ex = new String[names.size()];
    for (int i = 0; i < ex.length; i++) ex[i] = "ec.I0_" + names.get(i);
    ev.set("expr", ex);
    ev.set("data", "dset1");
    double[][] r = ev.getReal();   // [expr][time]
    double[] tt = m.sol("sol1").getPVals();
    double[] fac = new double[names.size()];
    for (int i = 0; i < fac.length; i++) fac[i] = names.get(i).endsWith("0") && names.get(i).length() == 2 ? 4 : 2;
    System.out.println("CSV,# full-detector channel currents in A for Q = 1e-15 C (quarter x 4 for X0/Y0, x 2 otherwise);"
                       + " Induce3D rho_s=" + rhoS + " case=" + kase + " t_ion_ns=" + tIon + " ve_um_ns=" + ve
                       + " order=" + order + " h_far=" + hFar + " h_near=" + hNear + " nz_gas=" + nzGas
                       + " elems=" + comp.mesh("mesh1").getNumElem());
    System.out.println("CSV,t_s," + String.join(",", names));
    for (int k = 0; k < tt.length; k++) {
      StringBuilder sb = new StringBuilder("CSV,").append(tt[k]);
      for (int i = 0; i < names.size(); i++) sb.append(',').append(fac[i] * r[i][k]);
      System.out.println(sb);
    }
    System.out.println("[induce3d] wrote " + tt.length + " times");
    return m;
  }

  /** pad (i, j) clipped to the quarter: {x0, x1, y0, y1} */
  static double[] pad(int i, int j, double P, double W) {
    return new double[]{Math.max(0, i * P - W / 2), i * P + W / 2, Math.max(0, j * P - W / 2), j * P + W / 2};
  }

  static void rect2(GeomSequence w, String tag, double x, double y, double dx, double dy) {
    GeomFeature f = w.create(tag, "Rectangle");
    f.set("pos", new double[]{x, y});
    f.set("size", new double[]{dx, dy});
  }

  static SelectionFeature box(ModelNode comp, String tag, int dim, double x0, double x1, double y0, double y1, double z0, double z1) {
    SelectionFeature s = comp.selection().create(tag, "Box");
    s.set("entitydim", dim);
    s.set("xmin", x0); s.set("xmax", x1); s.set("ymin", y0); s.set("ymax", y1); s.set("zmin", z0); s.set("zmax", z1);
    s.set("condition", "inside");
    return s;
  }

  static void union(ModelNode comp, String tag, int dim, List<String> in) {
    SelectionFeature u = comp.selection().create(tag, "Union");
    u.set("entitydim", dim);
    u.set("input", in.toArray(new String[0]));
  }

  static void diff(ModelNode comp, String tag, int dim, String add, String sub) {
    SelectionFeature d = comp.selection().create(tag, "Difference");
    d.set("entitydim", dim);
    d.set("add", new String[]{add});
    d.set("subtract", new String[]{sub});
  }

  static void sweepLayers(MeshSequence ms, String tag, String sel, int n) {
    MeshFeature d = ms.feature("swe").create(tag, "Distribution");
    d.selection().named(sel);
    d.set("numelem", n);
  }

  static void medium(Physics ec, String tag, String sel, String sigma, String epsr) {
    PhysicsFeature f = ec.create(tag, "CurrentConservation", 3);
    f.selection().named(sel);
    f.set("sigma_mat", "userdef");
    f.set("sigma", new String[]{sigma, "0", "0", "0", sigma, "0", "0", "0", sigma});
    f.set("epsilonr_mat", "userdef");
    f.set("epsilonr", new String[]{epsr, "0", "0", "0", epsr, "0", "0", "0", epsr});
  }

  public static void main(String[] args) throws Exception { run(); }
}
