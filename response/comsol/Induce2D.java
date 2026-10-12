/*
 * Induce2D.java -- direct time-domain induction in the MX17 readout stack (2D x-z cut).
 *
 * Independent check of the response chain's induction stage (S1/W2 kernels +
 * kernel_lut.apply_ion_transit). No weighting-field theorem is used: the avalanche
 * is a prescribed convection current Je = rho*v in the gas (electrons down onto the
 * ESL, ions up to the mesh), COMSOL solves charge conservation
 *     div( sigma E + eps dE/dt + Je ) = 0
 * through gas / resistive ESL / kapton / glue / pad copper / FR4, and the readout
 * pads are Terminals held at 0 V whose currents ARE the induced signals.
 *
 * 2D: the cut is across the ESL strips and pad columns (both run along y), so a
 * "pad" here is a continuous column, not the real checkerboard comb (that is 3D).
 * Geometry coordinates: x lateral, y vertical (COMSOL 2D), unit um, y = 0 at the
 * ESL top surface, mesh at y = GAP.
 *
 * Environment (all optional):
 *   RHO_S     sheet resistivity [ohm/sq]          default 2e6
 *   ESL_MODE  strips | uniform                     default strips
 *   X0_UM     avalanche position                   default 0
 *   CASE      pair | deposit                       default pair
 *             pair    = electron/ion pairs born on the exp(-z/L) profile at t = 0
 *             deposit = electrons born at the mesh, no ions: the chain's surface kernel G
 *   T_ION_NS  ion transit across the full gap      default 340
 *
 * Output: CSV lines on stdout prefixed "CSV,"; run_lxplus.sh writes them to $OUT.
 */
import com.comsol.model.*;
import com.comsol.model.util.*;
import com.comsol.model.physics.*;
import java.io.*;
import java.util.*;

public class Induce2D {
  static String env(String k, String d) { String v = System.getenv(k); return (v == null || v.isEmpty()) ? d : v; }

  public static Model run() throws Exception {
    double rhoS = Double.parseDouble(env("RHO_S", "2e6"));
    String eslMode = env("ESL_MODE", "strips");
    double x0 = Double.parseDouble(env("X0_UM", "0"));
    String kase = env("CASE", "pair");
    double tIon = Double.parseDouble(env("T_ION_NS", "340"));

    // stack [um] -- response/common/constants.py, solver/v6_pad_gaps.py
    final double GAP = 150, T_ESL = 10, T_KAP = 50, T_GLUE = 18.76, T_CU = 26, T_FR4 = 264;
    final double PAD_P = 780, PAD_W = 680, ESL_P = 800, ESL_W = 550, W = 6240; // +-8 pad pitches
    final int NPAD = 8, NTERM = 2;  // terminals: pads -2..2 around x = 0
    double yK = -T_ESL, yG = yK - T_KAP, yP = yG - T_GLUE, yT = yP - T_CU, yB = yT - T_FR4;

    Model m = ModelUtil.create("Model");
    m.param().set("rho_s", rhoS + "[ohm]");
    m.param().set("t_esl", T_ESL + "[um]");
    m.param().set("gap", GAP + "[um]");
    m.param().set("Lz", "14.8[um]", "ion-birth decay length 1/alpha (ions.py)");
    m.param().set("Lb", "5[um]", "deposit case: electrons born within this of the mesh");
    m.param().set("sig0", "33[um]", "avalanche lateral sigma (T7)");
    m.param().set("x0", x0 + "[um]");
    m.param().set("vi", "gap/" + tIon + "[ns]");
    // electron velocity slowed 10x from ~150 um/ns so the prompt pulse is resolvable; it only
    // stretches a ~1 ns spike to ~10 ns, invisible behind a 180 ns shaper
    m.param().set("ve", "15[um/ns]", "electron velocity (slowed 10x)");
    m.param().set("Q", "1e-12[C/m]", "line charge per unit depth; output is normalised by Q");
    m.param().set("ion_on", kase.equals("deposit") ? "0" : "1");
    m.param().set("dep", kase.equals("deposit") ? "1" : "0");

    ModelNode comp = m.component().create("comp1", true);
    GeomSequence g = comp.geom().create("geom1", 2);
    g.lengthUnit("um");
    g.selection().create("csel_esl", "CumulativeSelection");
    g.selection().create("csel_pad", "CumulativeSelection");
    rect(g, "gas", -W, 0, 2 * W, GAP, null);
    rect(g, "esllay", -W, yK, 2 * W, T_ESL, eslMode.equals("uniform") ? "csel_esl" : null);
    if (!eslMode.equals("uniform")) {
      int n = 0;
      for (double c = -Math.floor(W / ESL_P) * ESL_P; c <= W; c += ESL_P) {
        double a = Math.max(-W, c - ESL_W / 2), b = Math.min(W, c + ESL_W / 2);
        rect(g, "es" + (n++), a, yK, b - a, T_ESL, "csel_esl");
      }
    }
    rect(g, "kap", -W, yG, 2 * W, T_KAP, null);
    rect(g, "glue", -W, yP, 2 * W, T_GLUE, null);
    rect(g, "trench", -W, yT, 2 * W, T_CU, null);
    for (int i = -NPAD; i <= NPAD; i++) {
      double c = i * PAD_P, a = Math.max(-W, c - PAD_W / 2), b = Math.min(W, c + PAD_W / 2);
      rect(g, "pd" + (i + NPAD), a, yT, b - a, T_CU, "csel_pad");
    }
    rect(g, "fr4", -W, yB, 2 * W, T_FR4, null);
    g.run();

    // physics lives everywhere except inside the pad copper
    SelectionFeature all = comp.selection().create("s_all", "Explicit");
    all.all();
    SelectionFeature phys = comp.selection().create("s_phys", "Difference");
    phys.set("entitydim", 2);
    phys.set("add", new String[]{"s_all"});
    phys.set("subtract", new String[]{"geom1_csel_pad_dom"});
    SelectionFeature gasSel = box(comp, "s_gas", 2, -W - 1, W + 1, -0.01, GAP + 0.01);
    SelectionFeature eslGas = comp.selection().create("s_eslgap", "Difference");  // printed-layer gaps = gas
    eslGas.set("entitydim", 2);
    eslGas.set("add", new String[]{box(comp, "s_esllay", 2, -W - 1, W + 1, yK - 0.01, 0.01).tag()});
    eslGas.set("subtract", new String[]{"geom1_csel_esl_dom"});
    box(comp, "s_kap", 2, -W - 1, W + 1, yG - 0.01, yK + 0.01);
    box(comp, "s_glue", 2, -W - 1, W + 1, yP - 0.01, yG + 0.01);
    SelectionFeature tr = comp.selection().create("s_trench", "Difference");
    tr.set("entitydim", 2);
    tr.set("add", new String[]{box(comp, "s_trlay", 2, -W - 1, W + 1, yT - 0.01, yP + 0.01).tag()});
    tr.set("subtract", new String[]{"geom1_csel_pad_dom"});
    box(comp, "s_fr4", 2, -W - 1, W + 1, yB - 0.01, yT + 0.01);
    box(comp, "s_mesh", 1, -W - 1, W + 1, GAP - 0.01, GAP + 0.01);
    box(comp, "s_bot", 1, -W - 1, W + 1, yB - 0.01, yB + 0.01);
    List<String> groundPads = new ArrayList<>();
    for (int i = -NPAD; i <= NPAD; i++) {
      double c = i * PAD_P;
      box(comp, "s_pad" + (i + NPAD), 1, Math.max(-W, c - PAD_W / 2) - 0.01, Math.min(W, c + PAD_W / 2) + 0.01, yT - 0.01, yP + 0.01);
      if (Math.abs(i) > NTERM) groundPads.add("s_pad" + (i + NPAD));
    }
    SelectionFeature gnd = comp.selection().create("s_gnd", "Union");
    gnd.set("entitydim", 1);
    List<String> gl = new ArrayList<>(groundPads); gl.add("s_mesh"); gl.add("s_bot");
    gnd.set("input", gl.toArray(new String[0]));

    Physics ec = comp.physics().create("ec", "ConductiveMedia", "geom1");
    ec.selection().named("s_phys");
    ec.feature("cucns1").set("sigma_mat", "userdef");
    ec.feature("cucns1").set("sigma", new String[]{"0", "0", "0", "0", "0", "0", "0", "0", "0"});
    ec.feature("cucns1").set("epsilonr_mat", "userdef");
    ec.feature("cucns1").set("epsilonr", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    medium(ec, "cu_esl", "geom1_csel_esl_dom", "1/(rho_s*t_esl)", "3");
    medium(ec, "cu_kap", "s_kap", "0", "3.5");
    medium(ec, "cu_glue", "s_glue", "0", "3.2");
    medium(ec, "cu_tr", "s_trench", "0", "3.2");
    medium(ec, "cu_fr4", "s_fr4", "0", "4.4");
    ec.create("gnd1", "Ground", 1).selection().named("s_gnd");
    String[] tn = new String[2 * NTERM + 1];
    for (int i = -NTERM; i <= NTERM; i++) {
      String name = "p" + (i < 0 ? "m" + (-i) : String.valueOf(i));
      tn[i + NTERM] = name;
      PhysicsFeature t = ec.create("term_" + name, "Terminal", 1);
      t.selection().named("s_pad" + (i + NPAD));
      t.set("TerminalType", "Voltage");
      t.set("V0", "0");
      t.set("TerminalName", name);
    }
    // the avalanche as a convection current, vertical only (y up). The abs() keeps the
    // exponent bounded outside the window: 0*exp(+800) is NaN, not 0, and killed the solver.
    m.func().create("pz", "Analytic");
    m.func("pz").set("funcname", "pz");
    m.func("pz").set("expr", "(s>=0)*(s<=gap)*exp(-abs(s)/Lz)/(Lz*(1-exp(-gap/Lz)))");
    m.func("pz").set("args", new String[]{"s"});
    m.func("pz").set("argunit", "m");
    m.func("pz").set("fununit", "1/m");
    m.func().create("pb", "Analytic");
    m.func("pb").set("funcname", "pb");
    m.func("pb").set("expr", "(s>=0)*(s<=gap)*exp(-abs(gap-s)/Lb)/Lb");
    m.func("pb").set("args", new String[]{"s"});
    m.func("pb").set("argunit", "m");
    m.func("pb").set("fununit", "1/m");
    comp.variable().create("var1");
    comp.variable("var1").selection().geom("geom1", 2);
    comp.variable("var1").selection().named("s_gas");
    comp.variable("var1").set("gx", "exp(-(x-x0)^2/(2*sig0^2))/(sqrt(2*pi)*sig0)");
    comp.variable("var1").set("Jy", "Q*gx*(ion_on*vi*pz(y-vi*t) + ve*((1-dep)*pz(y+ve*t) + dep*pb(y+ve*t)))");
    PhysicsFeature je = ec.create("ecd1", "ExternalCurrentDensity", 2);
    je.selection().named("s_gas");
    je.set("Je", new String[]{"0", "Jy", "0"});

    // mesh: fine in the avalanche column and along the ESL, coarser elsewhere
    MeshSequence ms = comp.mesh().create("mesh1");
    box(comp, "s_col", 2, -W - 1, W + 1, -0.01, GAP + 0.01);
    ms.create("sz", "Size").set("hmax", "30[um]");
    ms.feature("sz").set("custom", "on");
    ms.feature("sz").set("hmaxactive", true);
    ms.feature("sz").set("hmin", "0.2[um]");
    ms.feature("sz").set("hminactive", true);
    ms.create("ftri", "FreeTri");
    ms.feature("ftri").create("size1", "Size");
    ms.feature("ftri").feature("size1").selection().geom("geom1", 2);
    ms.feature("ftri").feature("size1").selection().named("s_gas");
    ms.feature("ftri").feature("size1").set("custom", "on");
    ms.feature("ftri").feature("size1").set("hmax", "8[um]");
    ms.feature("ftri").feature("size1").set("hmaxactive", true);
    ms.run();
    System.out.println("[induce2d] mesh elements: " + comp.mesh("mesh1").getNumElem());

    // time-dependent solve, 0..3 us
    m.study().create("std1");
    m.study("std1").create("time", "Transient");
    m.study("std1").feature("time").set("tlist", "range(0,0.1[ns],2[ns]) range(3[ns],1[ns],3000[ns])");
    m.study("std1").createAutoSequences("all");
    SolverFeature td = m.sol("sol1").feature("t1");
    td.set("tstepsbdf", "free");
    td.set("maxstepconstraintbdf", "const");
    td.set("maxstepbdf", "1[ns]");
    td.set("initialstepbdfactive", true);
    td.set("initialstepbdf", "1e-13");
    m.sol("sol1").feature("t1").set("rtol", "1e-4");
    // V is ~10 mV for Q = 1 pC/m; automatic scaling from an all-zero initial state is degenerate
    m.sol("sol1").feature("v1").set("scalemethod", "manual");
    m.sol("sol1").feature("v1").set("scaleval", "1e-2");
    m.sol("sol1").runAll();

    // currents of every terminal, written as CSV (A per metre of depth, per C/m)
    NumericalFeature ev = m.result().numerical().create("gev1", "EvalGlobal");
    String[] ex = new String[tn.length];
    for (int i = 0; i < tn.length; i++) ex[i] = "ec.I0_" + tn[i];
    ev.set("expr", ex);
    ev.set("data", "dset1");
    double[][] r = ev.getReal();   // [expr][time]
    double[] tt = m.sol("sol1").getPVals();
    // COMSOL's file-system security preference refuses writes outside its own dirs, so the
    // CSV goes to stdout with a "CSV," prefix and run_lxplus.sh pulls it out of the log.
    System.out.println("CSV,# currents in A per metre of depth for Q = 1e-12 C/m; Induce2D rho_s=" + rhoS + " esl=" + eslMode + " x0_um=" + x0 + " case=" + kase + " t_ion_ns=" + tIon
                       + " elems=" + comp.mesh("mesh1").getNumElem());
    System.out.println("CSV,t_s," + String.join(",", tn));
    for (int k = 0; k < tt.length; k++) {
      StringBuilder sb = new StringBuilder("CSV,").append(tt[k]);
      for (int i = 0; i < tn.length; i++) sb.append(',').append(r[i][k]);
      System.out.println(sb);
    }
    System.out.println("[induce2d] wrote " + tt.length + " times");
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
