/*
 * InduceAxi.java -- direct time-domain induction, 2D-axisymmetric (r, z).
 *
 * The follow-up to Induce2D. In a 2D cut each ESL strip is an infinite floating line,
 * so charge that lands on it can neither spread nor drain and the surface kernel G is
 * flat after landing. The open question (T14: data want f_eff ~ 0.17-0.31 against
 * f_ion = 0.90) is about a POINT charge spreading on the resistive sheet during the
 * ~340 ns ion transit. The cheapest geometry that has that is axisymmetric: a point
 * avalanche on the axis above a uniform sheet, which spreads radially.
 *
 * Same method as Induce2D: the avalanche is a prescribed convection current
 * Je = rho*v (electrons down onto the ESL, ions up to the mesh), COMSOL solves
 *     div( sigma E + eps dE/dt + Je ) = 0
 * through gas / ESL / kapton / glue / pad copper / FR4, and the pads are Terminals at
 * 0 V whose currents are the induced signals. No weighting-field theorem is used.
 *
 * Pads: a central disk (r < PAD_W/2) and rings of width PAD_W on the PAD_P pitch,
 * standing in for the readout pads at their real scale (the real comb channels need 3D).
 * The sheet is grounded at its rim, R = 8 mm, far beyond the spreading reach in 3 us.
 *
 * Environment (all optional):
 *   RHO_S     sheet resistivity [ohm/sq]          default 2e6
 *   CASE      pair | deposit                       default pair
 *             pair    = electron/ion pairs born on the exp(-z/L) profile at t = 0
 *             deposit = electrons born at the mesh, no ions: the chain's surface kernel G
 *   T_ION_NS  ion transit across the full gap      default 340
 *   T_END_NS  end of the solve                     default 3000
 *   VE_UM_NS  electron velocity [um/ns]            default 15 (real ~150; see below)
 *
 * Output: CSV lines on stdout prefixed "CSV,"; run_lxplus.sh writes them to $OUT.
 */
import com.comsol.model.*;
import com.comsol.model.util.*;
import com.comsol.model.physics.*;
import java.util.*;

public class InduceAxi {
  static String env(String k, String d) { String v = System.getenv(k); return (v == null || v.isEmpty()) ? d : v; }

  public static Model run() throws Exception {
    double rhoS = Double.parseDouble(env("RHO_S", "2e6"));
    String kase = env("CASE", "pair");
    double tIon = Double.parseDouble(env("T_ION_NS", "340"));
    double tEnd = Double.parseDouble(env("T_END_NS", "3000"));
    double ve = Double.parseDouble(env("VE_UM_NS", "15"));

    // stack [um] -- response/common/constants.py, solver/v6_pad_gaps.py
    final double GAP = 150, T_ESL = 10, T_KAP = 50, T_GLUE = 18.76, T_CU = 26, T_FR4 = 264;
    final double PAD_P = 780, PAD_W = 680, R = 8000;
    final int NTERM = 3;    // terminals: disk + 2 rings; further rings grounded
    double yK = -T_ESL, yG = yK - T_KAP, yP = yG - T_GLUE, yT = yP - T_CU, yB = yT - T_FR4;

    Model m = ModelUtil.create("Model");
    m.param().set("rho_s", rhoS + "[ohm]");
    m.param().set("t_esl", T_ESL + "[um]");
    m.param().set("gap", GAP + "[um]");
    m.param().set("Lz", "14.8[um]", "ion-birth decay length 1/alpha (ions.py)");
    m.param().set("Lb", "10[um]", "deposit case: electrons born within this of the mesh");
    m.param().set("sig0", "33[um]", "avalanche lateral sigma (T7)");
    m.param().set("vi", "gap/" + tIon + "[ns]");
    // electron velocity: default slowed 10x from ~150 um/ns so the prompt pulse is resolvable.
    // Not harmless at low rho_s: the sheet relaxes during a 10 ns landing and erodes G's peak,
    // so G for the chain comparison should come from VE_UM_NS=150.
    m.param().set("ve", ve + "[um/ns]", "electron velocity");
    m.param().set("Q", "1e-15[C]", "avalanche charge; output is normalised by Q");
    m.param().set("ion_on", kase.equals("deposit") ? "0" : "1");
    m.param().set("dep", kase.equals("deposit") ? "1" : "0");

    ModelNode comp = m.component().create("comp1", true);
    GeomSequence g = comp.geom().create("geom1", 2);
    g.axisymmetric(true);
    g.lengthUnit("um");
    g.selection().create("csel_pad", "CumulativeSelection");
    rect(g, "gas", 0, 0, R, GAP, null);
    rect(g, "gascol", 0, 0, 300, GAP, null);       // fine-mesh column around the axis
    rect(g, "esl", 0, yK, R, T_ESL, null);
    rect(g, "kap", 0, yG, R, T_KAP, null);
    rect(g, "glue", 0, yP, R, T_GLUE, null);
    rect(g, "trench", 0, yT, R, T_CU, null);
    List<double[]> pads = new ArrayList<>();
    pads.add(new double[]{0, PAD_W / 2});
    for (int k = 1; k * PAD_P + PAD_W / 2 <= R; k++) pads.add(new double[]{k * PAD_P - PAD_W / 2, k * PAD_P + PAD_W / 2});
    for (int k = 0; k < pads.size(); k++) {
      double[] p = pads.get(k);
      rect(g, "pd" + k, p[0], yT, p[1] - p[0], T_CU, "csel_pad");
    }
    rect(g, "fr4", 0, yB, R, T_FR4, null);
    g.run();

    // physics lives everywhere except inside the pad copper
    SelectionFeature all = comp.selection().create("s_all", "Explicit");
    all.all();
    SelectionFeature phys = comp.selection().create("s_phys", "Difference");
    phys.set("entitydim", 2);
    phys.set("add", new String[]{"s_all"});
    phys.set("subtract", new String[]{"geom1_csel_pad_dom"});
    box(comp, "s_gas", 2, -1, R + 1, -0.01, GAP + 0.01);
    box(comp, "s_col", 2, -1, 300.01, -0.01, GAP + 0.01);
    box(comp, "s_esl", 2, -1, R + 1, yK - 0.01, 0.01);
    box(comp, "s_kap", 2, -1, R + 1, yG - 0.01, yK + 0.01);
    box(comp, "s_glue", 2, -1, R + 1, yP - 0.01, yG + 0.01);
    SelectionFeature tr = comp.selection().create("s_trench", "Difference");
    tr.set("entitydim", 2);
    tr.set("add", new String[]{box(comp, "s_trlay", 2, -1, R + 1, yT - 0.01, yP + 0.01).tag()});
    tr.set("subtract", new String[]{"geom1_csel_pad_dom"});
    box(comp, "s_fr4", 2, -1, R + 1, yB - 0.01, yT + 0.01);
    box(comp, "s_mesh", 1, -1, R + 1, GAP - 0.01, GAP + 0.01);
    box(comp, "s_bot", 1, -1, R + 1, yB - 0.01, yB + 0.01);
    box(comp, "s_rim", 1, R - 0.01, R + 0.01, yK - 0.01, 0.01);   // ESL rim to ground
    List<String> groundPads = new ArrayList<>();
    for (int k = 0; k < pads.size(); k++) {
      double[] p = pads.get(k);
      // the pad's boundary inside the physics: its top, sides and bottom
      box(comp, "s_pad" + k, 1, p[0] - 0.01, p[1] + 0.01, yT - 0.01, yP + 0.01);
      if (k >= NTERM) groundPads.add("s_pad" + k);
    }
    SelectionFeature gnd = comp.selection().create("s_gnd", "Union");
    gnd.set("entitydim", 1);
    List<String> gl = new ArrayList<>(groundPads); gl.add("s_mesh"); gl.add("s_bot"); gl.add("s_rim");
    gnd.set("input", gl.toArray(new String[0]));

    Physics ec = comp.physics().create("ec", "ConductiveMedia", "geom1");
    ec.selection().named("s_phys");
    ec.feature("cucns1").set("sigma_mat", "userdef");
    ec.feature("cucns1").set("sigma", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    ec.feature("cucns1").set("epsilonr_mat", "userdef");
    ec.feature("cucns1").set("epsilonr", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    medium(ec, "cu_esl", "s_esl", "1/(rho_s*t_esl)", "3");
    medium(ec, "cu_kap", "s_kap", "0", "3.5");
    medium(ec, "cu_glue", "s_glue", "0", "3.2");
    medium(ec, "cu_tr", "s_trench", "0", "3.2");
    medium(ec, "cu_fr4", "s_fr4", "0", "4.4");
    ec.create("gnd1", "Ground", 1).selection().named("s_gnd");
    String[] tn = new String[NTERM];
    for (int k = 0; k < NTERM; k++) {
      tn[k] = "p" + k;
      PhysicsFeature t = ec.create("term_" + tn[k], "Terminal", 1);
      t.selection().named("s_pad" + k);
      t.set("TerminalType", "Voltage");
      t.set("V0", "0");
      t.set("TerminalName", tn[k]);
    }
    // the avalanche as a convection current, vertical only (z up). The abs() keeps the
    // exponent bounded outside the window: 0*exp(+800) is NaN, not 0 (see Induce2D).
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
    comp.variable("var1").selection().geom("geom1", 2);
    comp.variable("var1").selection().named("s_gas");
    comp.variable("var1").set("gr", "exp(-r^2/(2*sig0^2))/(2*pi*sig0^2)");
    comp.variable("var1").set("Jz", "Q*gr*(ion_on*vi*pz(z-vi*t) + ve*((1-dep)*pz(z+ve*t) + dep*pb(z+ve*t)))");
    PhysicsFeature je = ec.create("ecd1", "ExternalCurrentDensity", 2);
    je.selection().named("s_gas");
    je.set("Je", new String[]{"0", "0", "Jz"});

    // mesh: fine in the avalanche column, graded outward
    MeshSequence ms = comp.mesh().create("mesh1");
    ms.create("sz", "Size");
    ms.feature("sz").set("custom", "on");
    ms.feature("sz").set("hmax", "40[um]");
    ms.feature("sz").set("hmaxactive", true);
    ms.feature("sz").set("hmin", "0.2[um]");
    ms.feature("sz").set("hminactive", true);
    ms.feature("sz").set("hgrad", 1.2);
    ms.feature("sz").set("hgradactive", true);
    ms.create("ftri", "FreeTri");
    ms.feature("ftri").create("size1", "Size");
    ms.feature("ftri").feature("size1").selection().geom("geom1", 2);
    ms.feature("ftri").feature("size1").selection().named("s_col");
    ms.feature("ftri").feature("size1").set("custom", "on");
    ms.feature("ftri").feature("size1").set("hmax", "2.5[um]");
    ms.feature("ftri").feature("size1").set("hmaxactive", true);
    ms.run();
    System.out.println("[induceaxi] mesh elements: " + comp.mesh("mesh1").getNumElem());

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
    // V is O(0.1 V) for Q = 1 fC; automatic scaling from an all-zero initial state is degenerate
    m.sol("sol1").feature("v1").set("scalemethod", "manual");
    m.sol("sol1").feature("v1").set("scaleval", "1e-1");
    m.sol("sol1").runAll();

    NumericalFeature ev = m.result().numerical().create("gev1", "EvalGlobal");
    String[] ex = new String[tn.length];
    for (int i = 0; i < tn.length; i++) ex[i] = "ec.I0_" + tn[i];
    ev.set("expr", ex);
    ev.set("data", "dset1");
    double[][] r = ev.getReal();   // [expr][time]
    double[] tt = m.sol("sol1").getPVals();
    // COMSOL's file-system security preference refuses writes outside its own dirs, so the
    // CSV goes to stdout with a "CSV," prefix and run_lxplus.sh pulls it out of the log.
    System.out.println("CSV,# currents in A for Q = 1e-15 C; InduceAxi rho_s=" + rhoS + " case=" + kase + " t_ion_ns=" + tIon + " ve_um_ns=" + ve
                       + " elems=" + comp.mesh("mesh1").getNumElem());
    System.out.println("CSV,t_s," + String.join(",", tn));
    for (int k = 0; k < tt.length; k++) {
      StringBuilder sb = new StringBuilder("CSV,").append(tt[k]);
      for (int i = 0; i < tn.length; i++) sb.append(',').append(r[i][k]);
      System.out.println(sb);
    }
    System.out.println("[induceaxi] wrote " + tt.length + " times");
    return m;
  }

  static void rect(GeomSequence g, String tag, double x, double y, double w, double h, String contrib) {
    GeomFeature f = g.create(tag, "Rectangle");
    f.set("pos", new double[]{x, y});
    f.set("size", new double[]{w, h});
    if (contrib != null) f.set("contributeto", contrib);
  }

  static SelectionFeature box(ModelNode comp, String tag, int dim, double x0, double x1, double y0, double y1) {
    SelectionFeature s = comp.selection().create(tag, "Box");
    s.set("entitydim", dim);
    s.set("xmin", x0); s.set("xmax", x1); s.set("ymin", y0); s.set("ymax", y1);
    s.set("condition", "inside");
    return s;
  }

  static void medium(Physics ec, String tag, String sel, String sigma, String epsr) {
    PhysicsFeature f = ec.create(tag, "CurrentConservation", 2);
    f.selection().named(sel);
    f.set("sigma_mat", "userdef");
    f.set("sigma", new String[]{sigma, "0", "0", "0", sigma, "0", "0", "0", sigma});
    f.set("epsilonr_mat", "userdef");
    f.set("epsilonr", new String[]{epsr, "0", "0", "0", epsr, "0", "0", "0", epsr});
  }

  public static void main(String[] args) throws Exception { run(); }
}
