import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import net.miginfocom.swing.MigLayout;

class GUI {
    public static JFrame frame;
    public static JCheckBox chkbx_showrc;
    public static JLabel lbl_stepone, lbl_steptwo, lbl_stepthree, lbl_stepfour, lbl_stepfive, lbl_stepsix,
            lbl_stepseven, lbl_r3phaseAct, lbl_r2phaseAct, lbl_r1phaseAct, lbl_r3phaseCalc, lbl_r2phaseCalc,
            lbl_r1phaseCalc, lbl_r3AphaseUnit, lbl_r2AphaseUnit, lbl_r1AphaseUnit, lbl_r3phaseUnit, lbl_r2phaseUnit,
            lbl_r1phaseUnit, lbl_estBkgPhaseVal, lbl_V0Val, lbl_d_V0Val, lbl_magMom, lbl_rho0, lbl_rho0val, lbl_dchi,
            lbl_dchiVal, lbl_a, lbl_aVal, lbl_MPcnt, lbl_Ri, lbl_ImRi, lbl_ReRi, lbl_echoDChi, lbl_aSE, lbl_err,
            lbl_errVal, lbl_gridSizeBase, lbl_gridSize, lbl_rho0SEVal, lbl_M, lbl_rcx, lbl_rcy, lbl_rcz,
            lbl_rczCorrection, lbl_eqPhase, lbl_eqPhaseUnit, lbl_rc, lbl_rcUnit, lbl_spx, lbl_spy, lbl_spz,
            lbl_spzCorrection, lbl_calculated, lbl_actual, lbl_r1, lbl_r1unit, lbl_r1phase, lbl_r2, lbl_r2phase, lbl_r3,
            lbl_r2unit, lbl_r3unit, lbl_r3phase, lbl_estBkgPhase, lbl_estBkgPhaseUnit, lbl_magMomUnit, lbl_snr,
            lbl_eps12, lbl_eps23, lbl_ReRiVal, lbl_ImRiVal, lbl_TEFirst, lbl_TEFirstUnit, lbl_RChi, lbl_RChiUnit,
            lbl_sigSE, lbl_spinCenter, lbl_innerBrack1, lbl_v1se, lbl_innerBrack2, lbl_v2se, lbl_outerBrack1,
            lbl_comma11, lbl_comma21, lbl_comma31, lbl_v2seYVal1, lbl_comma12, lbl_comma22, lbl_comma32,
            lbl_innerBrack3, lbl_outerBrack2, lbl_V0, lbl_rho0SE, lbl_V0Unit, lbl_echoDChiVal, lbl_aSEVal, lbl_comma23,
            lbl_comma33, lbl_comma24, lbl_comma34;
    public static JTextField txt_eqPhaseRC, txt_rcx, txt_rcy, txt_rcz, txt_rc, txt_r1,
            txt_r2, txt_r3, txt_spx, txt_spy, txt_spz, txt_snrVal, txt_eps12val, txt_eps23val, B0Text,
            txt_RChiVal, txt_magMomVal, txt_Ri, txt_M, txt_TEFirstVal, TE_lastText, txt_spinCenterXVal,
            txt_spinCenterYVal, txt_spinCenterZVal, txt_v1seXVal1, txt_v1seYVal1, txt_v1seZVal1, txt_v1seXVal2,
            txt_v1seYVal2, txt_v1seZVal2, txt_v2seXVal1, txt_v2seYVal1, txt_v2seZVal1, txt_v2seXVal2, txt_v2seYVal2,
            txt_v2seZVal2, txt_sigSEVal;
    public static JButton btn_loadImages, btn_estCR, btn_genSubpix, btn_estSubC, btn_verifyRadii, btn_removeBkg,
            btn_estBkgDens, btn_loadTE, btn_unk, btn_loadspinecho, btn_estRadSpinEcho, btn_redraw, btn_plotX, btn_plotY,
            btn_plotZ, btn_calcMagMom, btn_loadSimImg, btn_sumRi;

    private static final String RHO = "\u03C1";
    private static final String LABEL_0 = "\u2080";
    private static final String EPSILON = "\u03B5";
    private static final String CHI = "\u03C7";
    private static final String _DELTA = "\u0394";
    private static final String DELTA = "\u03B4";
    private static final String ITALICIZED_I = "\uD835\uDC8A";
    private static final String _SIGMA = "\u03C3";
    private static final String LABEL_0SE = "0,SE";
    private static final String LABEL_1SE = "1,SE";
    private static final String LABEL_2SE = "2,SE";
    private static final String PLUS_MINUS = "\u00B1";
    private static final String CUBED = "\u00B3";

    public GUI() {
        initialize();
        frame.setVisible(true);
    }

    /*
     * Built using Eclipse and MigLayout. Much easier to manage then before.
     * If this is to be edited I advise against changing code directly.
     * Copy and paste this function into Eclipse (or any other IDE) and
     * make sure you have some GUI plug-in configured so that you can view it
     * without compiling + running every time. Also make sure to have MigLayout
     * configured.
     */
    private void initialize() {
        try {
            frame = new JFrame("Calculate Magnetic Moment 3D");
            frame.setAlwaysOnTop(false);
            frame.setBounds(100, 100, 800, 700);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.getContentPane().setLayout(
                    new MigLayout("",
                            "[][][grow][45px:45.00px:45px][45px:45px:45px,grow][45px:45px:45px,grow][48.00,grow][]",
                            "[][][][][][][][][][][][][][][][][][][][][][][]"));
        } catch (Exception e) {
            Calculate_Magnetic_Moment_3D.logger.addError(e.toString());
        }

        lbl_stepone = new JLabel("1.");
        frame.getContentPane().add(lbl_stepone, "cell 0 0");

        btn_loadImages = new JButton("Load Magnitude and Phase Images");
        frame.getContentPane().add(btn_loadImages, "cell 1 0,growx");
        btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_steptwo = new JLabel("2.");
        frame.getContentPane().add(lbl_steptwo, "cell 0 1");

        btn_estCR = new JButton("Estimate Center/Radii");
        frame.getContentPane().add(btn_estCR, "cell 1 1,growx");
        btn_estCR.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_M = new JLabel("|M|:");
        frame.getContentPane().add(lbl_M, "flowx,cell 2 1,alignx left");

        lbl_rcx = new JLabel("x=");
        frame.getContentPane().add(lbl_rcx, "flowx,cell 3 1");

        lbl_rcy = new JLabel("y=");
        frame.getContentPane().add(lbl_rcy, "flowx,cell 4 1,alignx left");

        lbl_rcz = new JLabel("z=");
        frame.getContentPane().add(lbl_rcz, "flowx,cell 5 1,alignx left");

        lbl_rczCorrection = new JLabel("-1");
        frame.getContentPane().add(lbl_rczCorrection, "cell 6 1,alignx left");

        lbl_eqPhase = new JLabel("Equatorial Phase at RCenter=");
        frame.getContentPane().add(lbl_eqPhase, "flowx,cell 1 2");

        lbl_eqPhaseUnit = new JLabel("radian(s)");
        frame.getContentPane().add(lbl_eqPhaseUnit, "cell 2 2");

        lbl_rc = new JLabel("RCenter=");
        frame.getContentPane().add(lbl_rc, "cell 3 2,alignx trailing");

        txt_rc = new JTextField();
        txt_rc.setColumns(3);
        frame.getContentPane().add(txt_rc, "cell 4 2,growx");

        lbl_rcUnit = new JLabel("pixels");
        frame.getContentPane().add(lbl_rcUnit, "cell 5 2");

        lbl_stepthree = new JLabel("3.");
        frame.getContentPane().add(lbl_stepthree, "cell 0 3");

        btn_genSubpix = new JButton("Generate Subpixel Grid/Data");
        frame.getContentPane().add(btn_genSubpix, "cell 1 3,growx");
        btn_genSubpix.addActionListener(new Calculate_Magnetic_Moment_3D());

        btn_removeBkg = new JButton("Remove Bkg");
        frame.getContentPane().add(btn_removeBkg, "cell 2 3");
        btn_removeBkg.addActionListener(new Calculate_Magnetic_Moment_3D());

        chkbx_showrc = new JCheckBox("Show RCenter");
        chkbx_showrc.setVerticalAlignment(SwingConstants.TOP);
        frame.getContentPane().add(chkbx_showrc, "cell 3 3");

        lbl_gridSize = new JLabel("Grid Size:");
        frame.getContentPane().add(lbl_gridSize, "flowx,cell 7 3");

        lbl_stepfour = new JLabel("4.");
        frame.getContentPane().add(lbl_stepfour, "cell 0 4");

        btn_estSubC = new JButton("Estimate Subpixel Center");
        frame.getContentPane().add(btn_estSubC, "cell 1 4,growx");
        btn_estSubC.addActionListener(new Calculate_Magnetic_Moment_3D());

        txt_rcx = new JTextField();
        frame.getContentPane().add(txt_rcx, "cell 3 1");
        txt_rcx.setColumns(4);

        txt_rcy = new JTextField();
        frame.getContentPane().add(txt_rcy, "cell 4 1");
        txt_rcy.setColumns(4);

        txt_rcz = new JTextField();
        frame.getContentPane().add(txt_rcz, "cell 5 1");
        txt_rcz.setColumns(4);

        txt_M = new JTextField();
        txt_M.setText("50");
        frame.getContentPane().add(txt_M, "cell 2 1,alignx left");
        txt_M.setColumns(2);

        lbl_MPcnt = new JLabel("%");
        frame.getContentPane().add(lbl_MPcnt, "cell 2 1");

        txt_eqPhaseRC = new JTextField();
        txt_eqPhaseRC.setText("1.0");
        frame.getContentPane().add(txt_eqPhaseRC, "cell 1 2");
        txt_eqPhaseRC.setColumns(5);

        lbl_spx = new JLabel("x=");
        frame.getContentPane().add(lbl_spx, "flowx,cell 3 4");

        txt_spx = new JTextField();
        lbl_spx.setLabelFor(txt_spx);
        frame.getContentPane().add(txt_spx, "cell 3 4");
        txt_spx.setColumns(4);

        lbl_spy = new JLabel("y=");
        frame.getContentPane().add(lbl_spy, "flowx,cell 4 4");

        txt_spy = new JTextField();
        lbl_spy.setLabelFor(txt_spy);
        txt_spy.setColumns(4);
        frame.getContentPane().add(txt_spy, "cell 4 4");

        lbl_spz = new JLabel("z=");
        frame.getContentPane().add(lbl_spz, "flowx,cell 5 4");

        txt_spz = new JTextField();
        lbl_spz.setLabelFor(txt_spz);
        txt_spz.setColumns(4);
        frame.getContentPane().add(txt_spz, "cell 5 4");

        lbl_spzCorrection = new JLabel("-1");
        lbl_spzCorrection.setLabelFor(txt_spz);
        frame.getContentPane().add(lbl_spzCorrection, "cell 6 4");

        btn_redraw = new JButton("Redraw Center");
        frame.getContentPane().add(btn_redraw, "flowx,cell 1 5");
        btn_redraw.addActionListener(new Calculate_Magnetic_Moment_3D());

        btn_verifyRadii = new JButton("Verify Radii");
        frame.getContentPane().add(btn_verifyRadii, "cell 1 5");
        btn_verifyRadii.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_calculated = new JLabel("Calculated");
        frame.getContentPane().add(lbl_calculated, "cell 3 5");

        lbl_actual = new JLabel("Actual");
        frame.getContentPane().add(lbl_actual, "cell 5 5");

        lbl_r1 = new JLabel("R1=");
        frame.getContentPane().add(lbl_r1, "flowx,cell 1 6");

        txt_r1 = new JTextField();
        frame.getContentPane().add(txt_r1, "cell 1 6");
        txt_r1.setColumns(5);

        lbl_r1unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r1unit, "cell 1 6");

        lbl_r1phase = new JLabel("R1 Corresponding Phase:");
        frame.getContentPane().add(lbl_r1phase, "cell 2 6");

        lbl_r1phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseCalc, "cell 3 6");

        lbl_r1phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r1phaseUnit, "cell 4 6");

        lbl_r1phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseAct, "cell 5 6");

        lbl_r1AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r1AphaseUnit, "cell 6 6");

        btn_plotX = new JButton("Plot X Phase Profiles");
        frame.getContentPane().add(btn_plotX, "cell 7 6");
        btn_plotX.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_r2 = new JLabel("R2=");
        frame.getContentPane().add(lbl_r2, "flowx,cell 1 7");

        lbl_r2phase = new JLabel("R2 Corresponding Phase:");
        frame.getContentPane().add(lbl_r2phase, "cell 2 7");

        lbl_r2phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseCalc, "cell 3 7");

        lbl_r2phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2phaseUnit, "cell 4 7");

        lbl_r2phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseAct, "cell 5 7");

        lbl_r2AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2AphaseUnit, "cell 6 7");

        btn_plotY = new JButton("Plot Y Phase Profiles");
        frame.getContentPane().add(btn_plotY, "cell 7 7");
        btn_plotY.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_r3 = new JLabel("R3=");
        frame.getContentPane().add(lbl_r3, "flowx,cell 1 8");

        txt_r2 = new JTextField();
        txt_r2.setColumns(5);
        frame.getContentPane().add(txt_r2, "cell 1 7");

        txt_r3 = new JTextField();
        txt_r3.setColumns(5);
        frame.getContentPane().add(txt_r3, "cell 1 8");

        lbl_r2unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r2unit, "cell 1 7");

        lbl_r3unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r3unit, "cell 1 8");

        lbl_r3phase = new JLabel("R3 Corresponding Phase:");
        frame.getContentPane().add(lbl_r3phase, "cell 2 8");

        lbl_r3phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseCalc, "cell 3 8");

        lbl_r3phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3phaseUnit, "cell 4 8");

        lbl_r3phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseAct, "cell 5 8");

        lbl_r3AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3AphaseUnit, "cell 6 8");

        lbl_gridSizeBase = new JLabel("10" + CUBED);
        frame.getContentPane().add(lbl_gridSizeBase, "cell 7 3");

        btn_plotZ = new JButton("Plot Z Phase Profiles");
        frame.getContentPane().add(btn_plotZ, "cell 7 8");
        btn_plotZ.addActionListener(new Calculate_Magnetic_Moment_3D());

        btn_estBkgDens = new JButton("Estimate Bkg & " + RHO + "0");
        frame.getContentPane().add(btn_estBkgDens, "cell 1 9,growx");
        btn_estBkgDens.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_rho0 = new JLabel(RHO + "0 =");
        frame.getContentPane().add(lbl_rho0, "flowx,cell 2 9");

        lbl_rho0val = new JLabel("");
        frame.getContentPane().add(lbl_rho0val, "cell 2 9");

        lbl_estBkgPhase = new JLabel("Estimated Background Phase =");
        frame.getContentPane().add(lbl_estBkgPhase, "cell 3 9");

        lbl_estBkgPhaseVal = new JLabel("");
        frame.getContentPane().add(lbl_estBkgPhaseVal, "cell 6 9");

        lbl_estBkgPhaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_estBkgPhaseUnit, "cell 7 9");

        lbl_stepfive = new JLabel("5.");
        frame.getContentPane().add(lbl_stepfive, "cell 0 10");

        btn_calcMagMom = new JButton("Calculate Magnetic Moment");
        frame.getContentPane().add(btn_calcMagMom, "cell 1 10,growx");
        btn_calcMagMom.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_magMom = new JLabel("|p|=");
        frame.getContentPane().add(lbl_magMom, "flowx,cell 2 10");

        btn_loadSimImg = new JButton("Load Simulated Images");
        frame.getContentPane().add(btn_loadSimImg, "cell 1 11,growx");
        btn_loadSimImg.addActionListener(new Calculate_Magnetic_Moment_3D());

        txt_magMomVal = new JTextField();
        frame.getContentPane().add(txt_magMomVal, "cell 2 10");
        txt_magMomVal.setColumns(4);

        lbl_magMomUnit = new JLabel("radians*pixel" + CUBED);
        frame.getContentPane().add(lbl_magMomUnit, "cell 2 10");

        lbl_snr = new JLabel("SNR=");
        frame.getContentPane().add(lbl_snr, "flowx,cell 1 12");

        txt_snrVal = new JTextField();
        txt_snrVal.setText("1.0");
        frame.getContentPane().add(txt_snrVal, "cell 1 12");
        txt_snrVal.setColumns(3);

        lbl_eps12 = new JLabel(EPSILON + "12 =");
        frame.getContentPane().add(lbl_eps12, "flowx,cell 1 13");

        txt_eps12val = new JTextField();
        txt_eps12val.setText("0.0");
        txt_eps12val.setColumns(3);
        frame.getContentPane().add(txt_eps12val, "cell 1 13");

        lbl_ReRi = new JLabel("Real(S" + ITALICIZED_I + ") =");
        frame.getContentPane().add(lbl_ReRi, "flowx,cell 2 13");

        lbl_eps23 = new JLabel(EPSILON + "23 =");
        frame.getContentPane().add(lbl_eps23, "flowx,cell 1 14");

        txt_eps23val = new JTextField();
        txt_eps23val.setText("0.0");
        txt_eps23val.setColumns(3);
        frame.getContentPane().add(txt_eps23val, "cell 1 14");

        lbl_err = new JLabel("      " + DELTA + RHO + "/" + RHO + " =");
        lbl_err.setHorizontalAlignment(SwingConstants.RIGHT);
        frame.getContentPane().add(lbl_err, "cell 1 13");

        lbl_errVal = new JLabel("");
        frame.getContentPane().add(lbl_errVal, "cell 1 13");

        lbl_Ri = new JLabel("R" + ITALICIZED_I + " =");
        frame.getContentPane().add(lbl_Ri, "flowx,cell 2 12");

        txt_Ri = new JTextField();
        txt_Ri.setColumns(4);
        frame.getContentPane().add(txt_Ri, "cell 2 12");

        lbl_ImRi = new JLabel("Imag(S" + ITALICIZED_I + ") =");
        frame.getContentPane().add(lbl_ImRi, "flowx,cell 2 14");

        lbl_ReRiVal = new JLabel("");
        frame.getContentPane().add(lbl_ReRiVal, "cell 2 13");

        lbl_ImRiVal = new JLabel("");
        frame.getContentPane().add(lbl_ImRiVal, "cell 2 14");

        btn_sumRi = new JButton("Sum");
        frame.getContentPane().add(btn_sumRi, "cell 2 12");
        btn_sumRi.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_stepsix = new JLabel("6.");
        frame.getContentPane().add(lbl_stepsix, "cell 0 15");

        btn_loadTE = new JButton("Load First TE Images");
        frame.getContentPane().add(btn_loadTE, "cell 1 15,growx");
        btn_loadTE.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_TEFirst = new JLabel("TE_first =");
        frame.getContentPane().add(lbl_TEFirst, "flowx,cell 2 15");

        lbl_dchi = new JLabel(_DELTA + CHI + " =");
        frame.getContentPane().add(lbl_dchi, "flowx,cell 4 15");

        btn_unk = new JButton("TODO: Name");
        frame.getContentPane().add(btn_unk, "cell 1 16,growx");
        btn_unk.addActionListener(new Calculate_Magnetic_Moment_3D());

        txt_TEFirstVal = new JTextField();
        frame.getContentPane().add(txt_TEFirstVal, "cell 2 15");
        txt_TEFirstVal.setColumns(4);

        lbl_TEFirstUnit = new JLabel("ms");
        frame.getContentPane().add(lbl_TEFirstUnit, "cell 2 15");

        lbl_RChi = new JLabel("      R" + _DELTA + CHI + " =");
        frame.getContentPane().add(lbl_RChi, "flowx,cell 2 16");

        txt_RChiVal = new JTextField();
        txt_RChiVal.setText("0.0");
        frame.getContentPane().add(txt_RChiVal, "cell 2 16");
        txt_RChiVal.setColumns(4);

        lbl_RChiUnit = new JLabel("pixels");
        frame.getContentPane().add(lbl_RChiUnit, "cell 2 16");

        lbl_dchiVal = new JLabel("");
        frame.getContentPane().add(lbl_dchiVal, "cell 4 15");

        lbl_a = new JLabel("  a =");
        frame.getContentPane().add(lbl_a, "flowx,cell 4 16");

        lbl_aVal = new JLabel("");
        frame.getContentPane().add(lbl_aVal, "cell 4 16");

        lbl_stepseven = new JLabel("7.");
        frame.getContentPane().add(lbl_stepseven, "cell 0 17");

        btn_loadspinecho = new JButton("Load Spin Echo Images");
        frame.getContentPane().add(btn_loadspinecho, "cell 1 17,growx");
        btn_loadspinecho.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_sigSE = new JLabel("        " + _SIGMA + "SE =");
        frame.getContentPane().add(lbl_sigSE, "flowx,cell 2 17,alignx left");

        txt_sigSEVal = new JTextField();
        frame.getContentPane().add(txt_sigSEVal, "cell 2 17");
        txt_sigSEVal.setColumns(4);

        lbl_spinCenter = new JLabel("Second Image Center:  (");
        frame.getContentPane().add(lbl_spinCenter, "flowx,cell 1 18,alignx trailing");

        txt_spinCenterXVal = new JTextField();
        txt_spinCenterXVal.setColumns(3);
        frame.getContentPane().add(txt_spinCenterXVal, "flowx,cell 2 18,alignx center");

        lbl_innerBrack1 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_innerBrack1, "flowx,cell 3 18");

        lbl_v1se = new JLabel("V1,SE Box Coordinates:  (");
        frame.getContentPane().add(lbl_v1se, "flowx,cell 1 19,alignx trailing");

        txt_v1seXVal1 = new JTextField();
        txt_v1seXVal1.setColumns(3);
        frame.getContentPane().add(txt_v1seXVal1, "flowx,cell 2 19,alignx center");

        lbl_innerBrack2 = new JLabel("-1 )        (");
        frame.getContentPane().add(lbl_innerBrack2, "cell 3 19,alignx left");

        txt_v1seXVal2 = new JTextField();
        frame.getContentPane().add(txt_v1seXVal2, "flowx,cell 4 19,alignx center");
        txt_v1seXVal2.setColumns(3);

        txt_v1seYVal2 = new JTextField();
        frame.getContentPane().add(txt_v1seYVal2, "flowx,cell 5 19,alignx center");
        txt_v1seYVal2.setColumns(3);

        txt_v1seZVal2 = new JTextField();
        txt_v1seZVal2.setColumns(3);
        frame.getContentPane().add(txt_v1seZVal2, "cell 6 19,alignx center");

        lbl_outerBrack1 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_outerBrack1, "cell 7 19");

        lbl_v2se = new JLabel("V2,SE Box Coordinates:  (");
        frame.getContentPane().add(lbl_v2se, "flowx,cell 1 20,alignx trailing");

        txt_v2seXVal1 = new JTextField();
        txt_v2seXVal1.setColumns(3);
        frame.getContentPane().add(txt_v2seXVal1, "flowx,cell 2 20,alignx center");

        lbl_comma11 = new JLabel(",");
        frame.getContentPane().add(lbl_comma11, "cell 2 18");

        lbl_comma21 = new JLabel(",");
        frame.getContentPane().add(lbl_comma21, "cell 2 19");

        lbl_comma31 = new JLabel(",");
        frame.getContentPane().add(lbl_comma31, "cell 2 20");

        txt_spinCenterYVal = new JTextField();
        txt_spinCenterYVal.setColumns(3);
        frame.getContentPane().add(txt_spinCenterYVal, "cell 2 18");

        txt_v1seYVal1 = new JTextField();
        txt_v1seYVal1.setColumns(3);
        frame.getContentPane().add(txt_v1seYVal1, "cell 2 19,alignx center");

        txt_v2seYVal1 = new JTextField();
        txt_v2seYVal1.setColumns(3);
        frame.getContentPane().add(txt_v2seYVal1, "cell 2 20");

        lbl_comma12 = new JLabel(",");
        frame.getContentPane().add(lbl_comma12, "cell 2 18");

        lbl_comma22 = new JLabel(",");
        frame.getContentPane().add(lbl_comma22, "cell 2 19");

        lbl_comma32 = new JLabel(",");
        frame.getContentPane().add(lbl_comma32, "cell 2 20");

        txt_spinCenterZVal = new JTextField();
        frame.getContentPane().add(txt_spinCenterZVal, "cell 2 18");
        txt_spinCenterZVal.setColumns(3);

        txt_v1seZVal1 = new JTextField();
        txt_v1seZVal1.setColumns(3);
        frame.getContentPane().add(txt_v1seZVal1, "cell 2 19");

        txt_v2seZVal1 = new JTextField();
        txt_v2seZVal1.setColumns(3);
        frame.getContentPane().add(txt_v2seZVal1, "cell 2 20");

        lbl_innerBrack3 = new JLabel("-1 )        (");
        frame.getContentPane().add(lbl_innerBrack3, "cell 3 20,alignx left");

        txt_v2seXVal2 = new JTextField();
        txt_v2seXVal2.setColumns(3);
        frame.getContentPane().add(txt_v2seXVal2, "flowx,cell 4 20,alignx center");

        txt_v2seYVal2 = new JTextField();
        txt_v2seYVal2.setColumns(3);
        frame.getContentPane().add(txt_v2seYVal2, "flowx,cell 5 20,alignx center");

        txt_v2seZVal2 = new JTextField();
        txt_v2seZVal2.setColumns(3);
        frame.getContentPane().add(txt_v2seZVal2, "cell 6 20,alignx center");

        lbl_outerBrack2 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_outerBrack2, "cell 7 20");

        btn_estRadSpinEcho = new JButton("Estimate Object Radius From Spin Echo");
        frame.getContentPane().add(btn_estRadSpinEcho, "cell 1 21");
        btn_estRadSpinEcho.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_V0 = new JLabel("V0 =");
        frame.getContentPane().add(lbl_V0, "flowx,cell 2 21");

        lbl_rho0SE = new JLabel(RHO + "0,SE =");
        frame.getContentPane().add(lbl_rho0SE, "flowx,cell 3 21");

        lbl_echoDChi = new JLabel(_DELTA + CHI + " =");
        frame.getContentPane().add(lbl_echoDChi, "flowx,cell 2 22");

        lbl_V0Val = new JLabel("     ");
        frame.getContentPane().add(lbl_V0Val, "cell 2 21");

        lbl_V0Unit = new JLabel("pixels" + CUBED);
        frame.getContentPane().add(lbl_V0Unit, "cell 2 21");

        lbl_echoDChiVal = new JLabel("");
        frame.getContentPane().add(lbl_echoDChiVal, "cell 2 22");

        lbl_rho0SEVal = new JLabel("");
        frame.getContentPane().add(lbl_rho0SEVal, "cell 3 21");

        lbl_aSE = new JLabel("a =");
        frame.getContentPane().add(lbl_aSE, "flowx,cell 3 22");

        lbl_aSEVal = new JLabel("");
        frame.getContentPane().add(lbl_aSEVal, "cell 3 22");

        lbl_comma23 = new JLabel(",");
        frame.getContentPane().add(lbl_comma23, "cell 4 19");

        lbl_comma33 = new JLabel(",");
        frame.getContentPane().add(lbl_comma33, "cell 4 20");

        lbl_comma24 = new JLabel(",");
        frame.getContentPane().add(lbl_comma24, "cell 5 19");

        lbl_comma34 = new JLabel(",");
        frame.getContentPane().add(lbl_comma34, "cell 5 20");
    }

}
