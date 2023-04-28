import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import java.awt.event.ActionListener;
import java.nio.channels.CancelledKeyException;
import java.awt.event.ActionEvent;
import net.miginfocom.swing.MigLayout;

class GUI {
    public JFrame frame;
    public JCheckBox chkbx_showrc;
    public JLabel lbl_stepone, lbl_steptwo, lbl_stepthree, lbl_stepfour, lbl_stepfive, lbl_stepsix,
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
    public DefTextField txt_eqPhaseRC, txt_rcx, txt_rcy, txt_rcz, txt_rc, txt_r1,
            txt_r2, txt_r3, txt_spx, txt_spy, txt_spz, txt_snrVal, txt_eps12val, txt_eps23val, txt_B0Val,
            txt_RChiVal, txt_magMomVal, txt_Ri, txt_M, txt_TEFirstVal, txt_TELastVal, txt_spinCenterXVal,
            txt_spinCenterYVal, txt_spinCenterZVal, txt_v1seXVal1, txt_v1seYVal1, txt_v1seZVal1, txt_v1seXVal2,
            txt_v1seYVal2, txt_v1seZVal2, txt_v2seXVal1, txt_v2seYVal1, txt_v2seZVal1, txt_v2seXVal2, txt_v2seYVal2,
            txt_v2seZVal2, txt_sigSEVal;
    public JButton btn_loadImages, btn_estCR, btn_genSubpix, btn_estSubC, btn_verifyRadii, btn_removeBkg,
            btn_estBkgDens, btn_loadTE, btn_unk, btn_loadspinecho, btn_estRadSpinEcho, btn_redraw, btn_plotX, btn_plotY,
            btn_plotZ, btn_calcMagMom, btn_loadSimImg, btn_sumRi;

    private final String ITALICIZED_I = "\uD835\uDC8A";

    public GUI() {
        this.initialize();
        this.frame.setVisible(true);
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
        frame = new JFrame("Calculate Magnetic Moment 3D");
        frame.setAlwaysOnTop(false);
        frame.setBounds(100, 100, 900, 700);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setLayout(
                new MigLayout("",
                        "[][][grow][][45px:45.00px:45px][45px:45px:45px,grow][45px:45px:45px,grow][48.00,grow][]",
                        "[][][][][][][][][][][][][][][][][][][][][][][][]"));

        lbl_stepone = new JLabel("1.");
        frame.getContentPane().add(lbl_stepone, "cell 0 0");

        btn_loadImages = new JButton("Load Magnitude and Phase Images");
        frame.getContentPane().add(btn_loadImages, "cell 1 0,growx");
        // btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_steptwo = new JLabel("2.");
        frame.getContentPane().add(lbl_steptwo, "cell 0 1");

        btn_estCR = new JButton("Estimate Center/Radii");
        frame.getContentPane().add(btn_estCR, "cell 1 1,growx");

        lbl_M = new JLabel("|M%|:");
        frame.getContentPane().add(lbl_M, "flowx,cell 4 1,alignx right");

        txt_M = new DefTextField("50");
        // txt_M.setText("50");
        frame.getContentPane().add(txt_M, "cell 5 1,alignx left");
        txt_M.setColumns(2);

        lbl_eqPhase = new JLabel("Equatorial Phase at RCenter=");
        frame.getContentPane().add(lbl_eqPhase, "flowx,cell 1 2");

        lbl_eqPhaseUnit = new JLabel("radian(s)");
        frame.getContentPane().add(lbl_eqPhaseUnit, "cell 2 2");

        lbl_rc = new JLabel("RCenter=");
        frame.getContentPane().add(lbl_rc, "cell 4 2,alignx trailing");

        txt_rc = new DefTextField();
        txt_rc.setColumns(3);
        frame.getContentPane().add(txt_rc, "cell 5 2,growx");

        lbl_rcUnit = new JLabel("pixels");
        frame.getContentPane().add(lbl_rcUnit, "cell 6 2");

        lbl_stepthree = new JLabel("3.");
        frame.getContentPane().add(lbl_stepthree, "cell 0 3");

        btn_genSubpix = new JButton("Generate Subpixel Grid/Data");
        frame.getContentPane().add(btn_genSubpix, "cell 1 3,growx");

        btn_removeBkg = new JButton("Remove Bkg");
        frame.getContentPane().add(btn_removeBkg, "cell 2 3");

        chkbx_showrc = new JCheckBox("Show RCenter");
        chkbx_showrc.setVerticalAlignment(SwingConstants.TOP);
        frame.getContentPane().add(chkbx_showrc, "cell 4 3");

        lbl_gridSize = new JLabel("<html>Grid Size: 10<sup>3</sup></html>");
        frame.getContentPane().add(lbl_gridSize, "flowx,cell 8 3");

        lbl_stepfour = new JLabel("4.");
        frame.getContentPane().add(lbl_stepfour, "cell 0 4");

        btn_estSubC = new JButton("Estimate Subpixel Center");
        frame.getContentPane().add(btn_estSubC, "cell 1 4,growx");

        txt_eqPhaseRC = new DefTextField("1.0");
        frame.getContentPane().add(txt_eqPhaseRC, "cell 1 2");
        txt_eqPhaseRC.setColumns(5);

        lbl_spx = new JLabel("x=");
        frame.getContentPane().add(lbl_spx, "flowx,cell 2 4");

        btn_redraw = new JButton("Redraw Center");
        frame.getContentPane().add(btn_redraw, "flowx,cell 1 5");

        btn_verifyRadii = new JButton("Verify Radii");
        frame.getContentPane().add(btn_verifyRadii, "cell 1 5");

        lbl_calculated = new JLabel("Calculated");
        frame.getContentPane().add(lbl_calculated, "cell 4 5");

        lbl_actual = new JLabel("Actual");
        frame.getContentPane().add(lbl_actual, "cell 6 5");

        lbl_r1 = new JLabel("R1=");
        frame.getContentPane().add(lbl_r1, "flowx,cell 1 6");

        txt_r1 = new DefTextField();
        frame.getContentPane().add(txt_r1, "cell 1 6");
        txt_r1.setColumns(5);

        lbl_r1unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r1unit, "cell 1 6");

        lbl_r1phase = new JLabel("R1 Corresponding Phase:");
        frame.getContentPane().add(lbl_r1phase, "cell 2 6");

        lbl_r1phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseCalc, "cell 4 6");

        lbl_r1phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r1phaseUnit, "cell 5 6");

        lbl_r1phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseAct, "cell 6 6");

        lbl_r1AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r1AphaseUnit, "cell 7 6");

        btn_plotX = new JButton("Plot X Phase Profiles");
        frame.getContentPane().add(btn_plotX, "cell 8 6");

        lbl_r2 = new JLabel("R2=");
        frame.getContentPane().add(lbl_r2, "flowx,cell 1 7");

        lbl_r2phase = new JLabel("R2 Corresponding Phase:");
        frame.getContentPane().add(lbl_r2phase, "cell 2 7");

        lbl_r2phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseCalc, "cell 4 7");

        lbl_r2phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2phaseUnit, "cell 5 7");

        lbl_r2phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseAct, "cell 6 7");

        lbl_r2AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2AphaseUnit, "cell 7 7");

        btn_plotY = new JButton("Plot Y Phase Profiles");
        frame.getContentPane().add(btn_plotY, "cell 8 7");

        lbl_r3 = new JLabel("R3=");
        frame.getContentPane().add(lbl_r3, "flowx,cell 1 8");

        txt_r2 = new DefTextField();
        txt_r2.setColumns(5);
        frame.getContentPane().add(txt_r2, "cell 1 7");

        txt_r3 = new DefTextField();
        txt_r3.setColumns(5);
        frame.getContentPane().add(txt_r3, "cell 1 8");

        lbl_r2unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r2unit, "cell 1 7");

        lbl_r3unit = new JLabel("pixels");
        frame.getContentPane().add(lbl_r3unit, "cell 1 8");

        lbl_r3phase = new JLabel("R3 Corresponding Phase:");
        frame.getContentPane().add(lbl_r3phase, "cell 2 8");

        lbl_r3phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseCalc, "cell 4 8");

        lbl_r3phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3phaseUnit, "cell 5 8");

        lbl_r3phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseAct, "cell 6 8");

        lbl_r3AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3AphaseUnit, "cell 7 8");

        btn_plotZ = new JButton("Plot Z Phase Profiles");
        frame.getContentPane().add(btn_plotZ, "cell 8 8");

        btn_estBkgDens = new JButton("<html>Estimate Bkg & &rho;0");
        frame.getContentPane().add(btn_estBkgDens, "cell 1 9,growx");

        lbl_rho0 = new JLabel("<html>&rho;0 =</html>");
        frame.getContentPane().add(lbl_rho0, "flowx,cell 2 9");

        lbl_rho0val = new JLabel("");
        frame.getContentPane().add(lbl_rho0val, "cell 2 9");

        lbl_estBkgPhase = new JLabel("Estimated Background Phase =");
        frame.getContentPane().add(lbl_estBkgPhase, "cell 4 9");

        lbl_estBkgPhaseVal = new JLabel("");
        frame.getContentPane().add(lbl_estBkgPhaseVal, "cell 7 9");

        lbl_estBkgPhaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_estBkgPhaseUnit, "cell 8 9");

        lbl_stepfive = new JLabel("5.");
        frame.getContentPane().add(lbl_stepfive, "cell 0 10");

        btn_calcMagMom = new JButton("Calculate Magnetic Moment");
        frame.getContentPane().add(btn_calcMagMom, "cell 1 10,growx");

        lbl_magMom = new JLabel("|p|=");
        frame.getContentPane().add(lbl_magMom, "flowx,cell 2 10");

        btn_loadSimImg = new JButton("Load Simulated Images");
        frame.getContentPane().add(btn_loadSimImg, "cell 1 11,growx");

        txt_magMomVal = new DefTextField();
        frame.getContentPane().add(txt_magMomVal, "cell 2 10");
        txt_magMomVal.setColumns(4);

        lbl_magMomUnit = new JLabel("<html>radians*pixel<sup>3</sup></html>");
        frame.getContentPane().add(lbl_magMomUnit, "cell 2 10");

        lbl_snr = new JLabel("SNR=");
        frame.getContentPane().add(lbl_snr, "flowx,cell 1 12");

        txt_snrVal = new DefTextField("1.0");
        frame.getContentPane().add(txt_snrVal, "cell 1 12");
        txt_snrVal.setColumns(3);

        lbl_eps12 = new JLabel("<html>&epsilon;12 =</html>");
        frame.getContentPane().add(lbl_eps12, "flowx,cell 1 13");

        txt_eps12val = new DefTextField();
        txt_eps12val.setColumns(3);
        frame.getContentPane().add(txt_eps12val, "cell 1 13");

        lbl_ReRi = new JLabel("Real(S" + ITALICIZED_I + ") = ");
        frame.getContentPane().add(lbl_ReRi, "flowx,cell 2 13");

        lbl_eps23 = new JLabel("<html>&epsilon;23 =</html>");
        frame.getContentPane().add(lbl_eps23, "flowx,cell 1 14");

        txt_eps23val = new DefTextField();
        txt_eps23val.setColumns(3);
        frame.getContentPane().add(txt_eps23val, "cell 1 14");

        lbl_err = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&delta&rho/&rho =</html>");
        lbl_err.setHorizontalAlignment(SwingConstants.RIGHT);
        frame.getContentPane().add(lbl_err, "cell 1 13");

        lbl_errVal = new JLabel("");
        frame.getContentPane().add(lbl_errVal, "cell 1 13");

        lbl_Ri = new JLabel("R" + ITALICIZED_I + " =");
        frame.getContentPane().add(lbl_Ri, "flowx,cell 2 12");

        txt_Ri = new DefTextField();
        txt_Ri.setColumns(4);
        frame.getContentPane().add(txt_Ri, "cell 2 12");

        lbl_ImRi = new JLabel("Imag(S" + ITALICIZED_I + ") = ");
        frame.getContentPane().add(lbl_ImRi, "flowx,cell 2 14");

        lbl_ReRiVal = new JLabel("");
        frame.getContentPane().add(lbl_ReRiVal, "cell 2 13");

        lbl_ImRiVal = new JLabel("");
        frame.getContentPane().add(lbl_ImRiVal, "cell 2 14");

        btn_sumRi = new JButton("Sum");
        frame.getContentPane().add(btn_sumRi, "cell 2 12");

        JLabel lbl_B0 = new JLabel("B0  =");
        frame.getContentPane().add(lbl_B0, "flowx,cell 1 15");

        txt_B0Val = new DefTextField();
        frame.getContentPane().add(txt_B0Val, "cell 1 15");
        txt_B0Val.setColumns(3);

        JLabel lbl_B0Unit = new JLabel("T");
        frame.getContentPane().add(lbl_B0Unit, "cell 1 15");

        JLabel lbl_TElast = new JLabel("    TE_last =");
        frame.getContentPane().add(lbl_TElast, "cell 1 15,alignx right");

        lbl_stepsix = new JLabel("6.");
        frame.getContentPane().add(lbl_stepsix, "cell 0 16");

        btn_loadTE = new JButton("Load First TE Images");
        frame.getContentPane().add(btn_loadTE, "cell 1 16,growx");

        lbl_TEFirst = new JLabel("TE_first =");
        frame.getContentPane().add(lbl_TEFirst, "flowx,cell 2 16");

        lbl_dchi = new JLabel("<html>&Delta;&Chi =</html>");
        frame.getContentPane().add(lbl_dchi, "flowx,cell 5 16");

        btn_unk = new JButton("TODO: Name");
        frame.getContentPane().add(btn_unk, "cell 1 17,growx");

        txt_TEFirstVal = new DefTextField();
        frame.getContentPane().add(txt_TEFirstVal, "cell 2 16");
        txt_TEFirstVal.setColumns(4);

        lbl_TEFirstUnit = new JLabel("ms");
        frame.getContentPane().add(lbl_TEFirstUnit, "cell 2 16");

        lbl_RChi = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;R<sub>&Delta&Chi</sub> =</html>");
        frame.getContentPane().add(lbl_RChi, "flowx,cell 2 17");

        txt_RChiVal = new DefTextField();
        frame.getContentPane().add(txt_RChiVal, "cell 2 17");
        txt_RChiVal.setColumns(4);

        lbl_RChiUnit = new JLabel("pixels");
        frame.getContentPane().add(lbl_RChiUnit, "cell 2 17");

        lbl_dchiVal = new JLabel("");
        frame.getContentPane().add(lbl_dchiVal, "cell 5 16");

        lbl_a = new JLabel("  a =");
        frame.getContentPane().add(lbl_a, "flowx,cell 5 17");

        lbl_aVal = new JLabel("");
        frame.getContentPane().add(lbl_aVal, "cell 5 17");

        lbl_stepseven = new JLabel("7.");
        frame.getContentPane().add(lbl_stepseven, "cell 0 18");

        btn_loadspinecho = new JButton("Load Spin Echo Images");
        frame.getContentPane().add(btn_loadspinecho, "cell 1 18,growx");

        lbl_sigSE = new JLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&sigma<sub>SE</sub> =</html>");
        frame.getContentPane().add(lbl_sigSE, "flowx,cell 2 18,alignx left");

        txt_sigSEVal = new DefTextField();
        frame.getContentPane().add(txt_sigSEVal, "cell 2 18");
        txt_sigSEVal.setColumns(4);

        lbl_spinCenter = new JLabel("Second Image Center:  (");
        frame.getContentPane().add(lbl_spinCenter, "flowx,cell 1 19,alignx trailing");

        txt_spinCenterXVal = new DefTextField();
        txt_spinCenterXVal.setColumns(4);
        frame.getContentPane().add(txt_spinCenterXVal, "flowx,cell 2 19,alignx left");

        // lbl_v1se = new JLabel("V1,SE Box Coordinates: (");
        lbl_v1se = new JLabel("<html>V<sub>1,SE</sub> Box Coordinates:  (</html>");
        frame.getContentPane().add(lbl_v1se, "flowx,cell 1 20,alignx trailing");

        txt_v1seXVal1 = new DefTextField();
        txt_v1seXVal1.setColumns(4);
        frame.getContentPane().add(txt_v1seXVal1, "flowx,cell 2 20,alignx left");

        JLabel lbl_brack23 = new JLabel("(");
        frame.getContentPane().add(lbl_brack23, "cell 3 20,alignx trailing");

        txt_v1seXVal2 = new DefTextField();
        frame.getContentPane().add(txt_v1seXVal2, "cell 4 20,alignx center");
        txt_v1seXVal2.setColumns(4);

        JLabel lbl_brack24 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_brack24, "cell 7 20");

        lbl_v2se = new JLabel("<html>V<sub>2,SE</sub> Box Coordinates:  (</html>");
        frame.getContentPane().add(lbl_v2se, "flowx,cell 1 21,alignx trailing");

        txt_v2seXVal1 = new DefTextField();
        txt_v2seXVal1.setColumns(4);
        frame.getContentPane().add(txt_v2seXVal1, "flowx,cell 2 21,alignx left");

        lbl_comma11 = new JLabel(",");
        frame.getContentPane().add(lbl_comma11, "cell 2 19");

        lbl_comma21 = new JLabel(",");
        frame.getContentPane().add(lbl_comma21, "cell 2 20");

        lbl_comma31 = new JLabel(",");
        frame.getContentPane().add(lbl_comma31, "cell 2 21");

        txt_spinCenterYVal = new DefTextField();
        txt_spinCenterYVal.setColumns(4);
        frame.getContentPane().add(txt_spinCenterYVal, "cell 2 19");

        txt_v1seYVal1 = new DefTextField();
        txt_v1seYVal1.setColumns(4);
        frame.getContentPane().add(txt_v1seYVal1, "cell 2 20,alignx center");

        txt_v2seYVal1 = new DefTextField();
        txt_v2seYVal1.setColumns(4);
        frame.getContentPane().add(txt_v2seYVal1, "cell 2 21");

        lbl_comma12 = new JLabel(",");
        frame.getContentPane().add(lbl_comma12, "cell 2 19");

        lbl_comma22 = new JLabel(",");
        frame.getContentPane().add(lbl_comma22, "cell 2 20");

        lbl_comma32 = new JLabel(",");
        frame.getContentPane().add(lbl_comma32, "cell 2 21");

        txt_spinCenterZVal = new DefTextField();
        frame.getContentPane().add(txt_spinCenterZVal, "cell 2 19");
        txt_spinCenterZVal.setColumns(4);

        txt_v1seZVal1 = new DefTextField();
        txt_v1seZVal1.setColumns(4);
        frame.getContentPane().add(txt_v1seZVal1, "cell 2 20");

        txt_v2seZVal1 = new DefTextField();
        txt_v2seZVal1.setColumns(4);
        frame.getContentPane().add(txt_v2seZVal1, "cell 2 21");

        JLabel lbl_brack33 = new JLabel("(");
        frame.getContentPane().add(lbl_brack33, "cell 3 21");

        txt_v2seXVal2 = new DefTextField();
        txt_v2seXVal2.setColumns(4);
        frame.getContentPane().add(txt_v2seXVal2, "cell 4 21,alignx center");

        JLabel lbl_brack34 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_brack34, "cell 7 21");

        btn_estRadSpinEcho = new JButton("Estimate Object Radius From Spin Echo");
        frame.getContentPane().add(btn_estRadSpinEcho, "cell 1 22");

        lbl_V0 = new JLabel("<html>V<sub>0</sub> =</html>");
        frame.getContentPane().add(lbl_V0, "flowx,cell 2 22");

        lbl_rho0SE = new JLabel("<html>&rho<sub>0,SE</sub> =</html>");
        frame.getContentPane().add(lbl_rho0SE, "flowx,cell 4 22");

        lbl_echoDChi = new JLabel("<html>&Delta&Chi =</html>");
        frame.getContentPane().add(lbl_echoDChi, "flowx,cell 2 23");

        lbl_V0Val = new JLabel("     ");
        frame.getContentPane().add(lbl_V0Val, "cell 2 22");

        lbl_V0Unit = new JLabel("<html>pixels<sup>3</sup></html>");
        frame.getContentPane().add(lbl_V0Unit, "cell 2 22,aligny top");

        lbl_echoDChiVal = new JLabel("");
        frame.getContentPane().add(lbl_echoDChiVal, "cell 2 23");

        lbl_rho0SEVal = new JLabel("");
        frame.getContentPane().add(lbl_rho0SEVal, "cell 4 22");

        lbl_aSE = new JLabel("a =");
        frame.getContentPane().add(lbl_aSE, "flowx,cell 4 23");

        lbl_aSEVal = new JLabel("");
        frame.getContentPane().add(lbl_aSEVal, "cell 4 23");

        lbl_comma23 = new JLabel(",");
        frame.getContentPane().add(lbl_comma23, "flowx,cell 5 20");

        lbl_comma33 = new JLabel(",");
        frame.getContentPane().add(lbl_comma33, "flowx,cell 5 21");

        lbl_comma24 = new JLabel(",");
        frame.getContentPane().add(lbl_comma24, "flowx,cell 6 20");

        lbl_comma34 = new JLabel(",");
        frame.getContentPane().add(lbl_comma34, "flowx,cell 6 21");
        lbl_spx.setLabelFor(txt_spx);

        txt_spx = new DefTextField();
        frame.getContentPane().add(txt_spx, "cell 2 4");
        txt_spx.setColumns(4);

        lbl_spy = new JLabel("y=");
        frame.getContentPane().add(lbl_spy, "cell 2 4");

        txt_spy = new DefTextField();
        txt_spy.setColumns(4);
        frame.getContentPane().add(txt_spy, "cell 2 4");
        lbl_spy.setLabelFor(txt_spy);

        lbl_spz = new JLabel("z=");
        frame.getContentPane().add(lbl_spz, "cell 2 4");

        txt_spz = new DefTextField();
        txt_spz.setColumns(4);
        frame.getContentPane().add(txt_spz, "cell 2 4");
        lbl_spz.setLabelFor(txt_spz);

        lbl_spzCorrection = new JLabel("-1");
        frame.getContentPane().add(lbl_spzCorrection, "cell 2 4");
        lbl_spzCorrection.setLabelFor(txt_spz);

        lbl_innerBrack1 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_innerBrack1, "cell 2 19");

        JLabel lbl_brack22 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_brack22, "cell 2 20,alignx left");

        JLabel lbl_brack32 = new JLabel("-1 )");
        frame.getContentPane().add(lbl_brack32, "cell 2 21,alignx left");

        txt_v1seYVal2 = new DefTextField();
        frame.getContentPane().add(txt_v1seYVal2, "cell 5 20,alignx center");
        txt_v1seYVal2.setColumns(4);

        txt_v2seYVal2 = new DefTextField();
        txt_v2seYVal2.setColumns(4);
        frame.getContentPane().add(txt_v2seYVal2, "cell 5 21,alignx center");

        txt_v1seZVal2 = new DefTextField();
        txt_v1seZVal2.setColumns(4);
        frame.getContentPane().add(txt_v1seZVal2, "cell 6 20,alignx center");

        txt_v2seZVal2 = new DefTextField();
        txt_v2seZVal2.setColumns(4);
        frame.getContentPane().add(txt_v2seZVal2, "cell 6 21,alignx center");

        lbl_rcx = new JLabel("x=");
        frame.getContentPane().add(lbl_rcx, "flowx,cell 2 1");

        txt_rcx = new DefTextField();
        frame.getContentPane().add(txt_rcx, "cell 2 1");
        txt_rcx.setColumns(4);

        lbl_rcy = new JLabel("y=");
        frame.getContentPane().add(lbl_rcy, "cell 2 1,alignx left");

        txt_rcy = new DefTextField();
        frame.getContentPane().add(txt_rcy, "cell 2 1");
        txt_rcy.setColumns(4);

        lbl_rcz = new JLabel("z=");
        frame.getContentPane().add(lbl_rcz, "cell 2 1,alignx left");

        txt_rcz = new DefTextField();
        frame.getContentPane().add(txt_rcz, "cell 2 1");
        txt_rcz.setColumns(4);

        lbl_rczCorrection = new JLabel("-1");
        frame.getContentPane().add(lbl_rczCorrection, "cell 2 1,alignx left");

        txt_TELastVal = new DefTextField();
        frame.getContentPane().add(txt_TELastVal, "cell 1 15,alignx right");
        txt_TELastVal.setColumns(3);

        JLabel lbl_TELastUnit = new JLabel("ms");
        lbl_TELastUnit.setHorizontalAlignment(SwingConstants.RIGHT);
        frame.getContentPane().add(lbl_TELastUnit, "cell 1 15,alignx right");
        frame.getContentPane().add(lbl_TELastUnit, "cell 1 15,alignx right");

        // ActionListeners
        btn_loadImages.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.load_mag_phase_images();
            }
        });
        btn_estCR.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.est_center_rad();
            }
        });
        btn_genSubpix.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.gen_subpix();
            }
        });
        btn_removeBkg.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.remove_bkg();
            }
        });
        btn_estSubC.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.est_subpix_ctr();
            }
        });
        btn_redraw.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.redraw_center();
            }
        });
        btn_verifyRadii.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.verify_radii();
            }
        });
        btn_plotX.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plot_x();
            }
        });
        btn_plotY.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plot_y();
            }
        });
        btn_plotZ.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plot_z();
            }
        });
        btn_estBkgDens.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.est_bkg_dens();
            }
        });
        btn_calcMagMom.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.calc_mag_moment();
            }
        });
        btn_loadSimImg.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.load_sim_images();
            }
        });
        btn_sumRi.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.sum_ri();
            }
        });
        btn_loadTE.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.load_first_TEImgs();
            }
        });
        btn_unk.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.noname();
            }
        });
        btn_loadspinecho.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.load_spin_echoImgs();
            }
        });
        btn_estRadSpinEcho.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.est_radius_spin_echo();
            }
        });
    }
}
