import java.awt.event.ActionListener;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.JTextField;
import net.miginfocom.swing.MigLayout;
import javax.swing.event.DocumentListener;

public class GUI {

    public JFrame frame;
    public JLabel lbl_stepone, lbl_steptwo, lbl_stepthree, lbl_stepfour, lbl_stepfive, lbl_stepsix,
            lbl_stepseven, lbl_r3phaseAct, lbl_r2phaseAct, lbl_r1phaseAct, lbl_r3phaseCalc, lbl_r2phaseCalc,
            lbl_r1phaseCalc, lbl_r3AphaseUnit, lbl_r2AphaseUnit, lbl_r1AphaseUnit, lbl_r3phaseUnit,
            lbl_r2phaseUnit,
            lbl_r1phaseUnit,
            lbl_calculated, lbl_actual, lbl_r1phase, lbl_r2phase,
            lbl_r3phase;
    public JButton btn_loadImages, btn_estCR, btn_genSubpix, btn_estSubC, btn_verifyRadii, btn_removeBkg,
            btn_estBkgDens, btn_showRCenter, btn_loadTE, btn_unk, btn_loadspinecho, btn_estRadSpinEcho, btn_redraw,
            btn_plotX, btn_plotY, btn_plotZ, btn_calcMagMom, btn_loadSimImg, btn_sumRi, btn_removeBkgS5;
    public LabeledTextField ltf_eqPhase, ltf_rc, ltf_r1, ltf_r2, ltf_r3, ltf_rcx, ltf_rcy, ltf_rcz, ltf_M, ltf_spx,
            ltf_spy, ltf_spz,
            ltf_TEFirst, ltf_TELast, ltf_B0, ltf_RChi, ltf_sigSE, ltf_eps12, ltf_eps23, ltf_snr, ltf_Ri, ltf_secondImgX,
            ltf_secondImgY, ltf_secondImgZ, ltf_v1seX1, ltf_v1seX2, ltf_v1seY1, ltf_v1seY2, ltf_v1seZ1, ltf_v1seZ2,
            ltf_v2seX1, ltf_v2seX2, ltf_v2seY1, ltf_v2seY2, ltf_v2seZ1, ltf_v2seZ2;
    public LabeledLabel ll_estBkgPhase, ll_grid, ll_rho0, ll_ReRi, ll_ImRi, ll_momenterror, ll_rho0SE, ll_aSE, ll_dChi,
            ll_dChiSE, ll_a, ll_V0, ltf_magMom;
    public LabeledDropDown<ImageItem.Axis> ldd_MriAxis;

    public GUI() {
        this.initialize();
        this.frame.setVisible(true);
    }

    /*
     * Built using Eclipse and MigLayout. Much easier to manage then before.
     * If this is to be edited I advise against changing code directly.
     * Copy and paste this file into Eclipse (or any other IDE) and
     * make sure you have some GUI plug-in configured so that you can view it
     * without compiling + running every time. Also make sure to have MigLayout
     * configured.
     */
    private void initialize() {
        frame = new JFrame("Calculate Magnetic Moment 3D");
        frame.setAlwaysOnTop(false);
        frame.setBounds(100, 100, 900, 850);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.getContentPane().setLayout(
                new MigLayout("",
                        "[][][grow][][45px:45.00px:45px][45px:45px:45px,grow][45px:64.00px:45px,grow][35.00,grow][]",
                        "[][][][][][][][][][][][][][][][][][][][][][][][][]"));
        Font font = new Font("Dialog", Font.PLAIN, 12);
        frame.setFont(font);

        lbl_stepone = new JLabel("1.");
        frame.getContentPane().add(lbl_stepone, "cell 0 0");

        btn_loadImages = new JButton("Load Magnitude and Phase Images");
        frame.getContentPane().add(btn_loadImages, "cell 1 0,growx");
        // btn_loadImages.addActionListener(new Calculate_Magnetic_Moment_3D());

        lbl_steptwo = new JLabel("2.");
        frame.getContentPane().add(lbl_steptwo, "cell 0 1");

        btn_estCR = new JButton("Estimate Center/Radii");
        frame.getContentPane().add(btn_estCR, "cell 1 1,growx");

        ltf_eqPhase = new LabeledTextField("Equatorial Phase at RCenter=", "1.0", "radian(s)", 5);
        frame.getContentPane().add(ltf_eqPhase, "cell 1 2");

        ltf_rc = new LabeledTextField("RCenter=", null, "pixels", 3);
        frame.getContentPane().add(ltf_rc, "cell 4 2 3 1");

        lbl_stepthree = new JLabel("3.");
        frame.getContentPane().add(lbl_stepthree, "cell 0 3");

        btn_genSubpix = new JButton("Generate Subpixel Grid/Data");
        frame.getContentPane().add(btn_genSubpix, "cell 1 3,growx");

        btn_showRCenter = new JButton("Show RCenter");
        frame.getContentPane().add(btn_showRCenter, "cell 4 3");

        ll_grid = new LabeledLabel("Grid Size:", "<html>10<sup>3</sup></html>", null);
        frame.getContentPane().add(ll_grid, "flowx, cell 8 3");

        lbl_stepfour = new JLabel("4.");
        frame.getContentPane().add(lbl_stepfour, "cell 0 4");

        btn_estSubC = new JButton("Estimate Subpixel Center");
        frame.getContentPane().add(btn_estSubC, "cell 1 4,growx");

        ltf_spx = new LabeledTextField("x=", "0.0", null, 4);
        frame.getContentPane().add(ltf_spx, "cell 2 4");

        ltf_spy = new LabeledTextField("y=", "0.0", null, 4);
        frame.getContentPane().add(ltf_spy, "cell 2 4");

        ltf_spz = new LabeledTextField("z=", "0.0", "-1", 4);
        frame.getContentPane().add(ltf_spz, "cell 2 4");

        btn_redraw = new JButton("Redraw Center");
        frame.getContentPane().add(btn_redraw, "flowx,cell 1 5");

        btn_verifyRadii = new JButton("Verify Radii");
        frame.getContentPane().add(btn_verifyRadii, "cell 1 5");

        lbl_calculated = new JLabel("Calculated");
        frame.getContentPane().add(lbl_calculated, "cell 4 5");

        lbl_actual = new JLabel("Actual");
        frame.getContentPane().add(lbl_actual, "cell 6 5,alignx right");

        ltf_r1 = new LabeledTextField("R1=", "0.0", "pixels", 4);
        frame.getContentPane().add(ltf_r1, "cell 1 6");

        lbl_r1phase = new JLabel("R1 Corresponding Phase:");
        frame.getContentPane().add(lbl_r1phase, "cell 2 6,alignx right");

        lbl_r1phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseCalc, "cell 4 6,alignx center");

        lbl_r1phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r1phaseUnit, "cell 5 6,alignx left");

        lbl_r1phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r1phaseAct, "cell 6 6,alignx right");

        lbl_r1AphaseUnit = new JLabel("radians");
        lbl_r1AphaseUnit.setHorizontalAlignment(SwingConstants.RIGHT);
        frame.getContentPane().add(lbl_r1AphaseUnit, "cell 7 6,alignx right");

        btn_plotX = new JButton("Plot X Phase Profiles");
        frame.getContentPane().add(btn_plotX, "cell 8 6");

        ltf_r2 = new LabeledTextField("R2=", "0.0", "pixels", 4);
        frame.getContentPane().add(ltf_r2, "cell 1 7");

        lbl_r2phase = new JLabel("R2 Corresponding Phase:");
        frame.getContentPane().add(lbl_r2phase, "cell 2 7,alignx right");

        lbl_r2phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseCalc, "cell 4 7,alignx center");

        lbl_r2phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2phaseUnit, "cell 5 7,alignx left");

        lbl_r2phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r2phaseAct, "cell 6 7,alignx right");

        lbl_r2AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r2AphaseUnit, "cell 7 7,alignx right");

        btn_plotY = new JButton("Plot Y Phase Profiles");
        frame.getContentPane().add(btn_plotY, "cell 8 7");

        ltf_r3 = new LabeledTextField("R3=", "0.0", "pixels", 4);
        frame.getContentPane().add(ltf_r3, "cell 1 8");

        lbl_r3phase = new JLabel("R3 Corresponding Phase:");
        frame.getContentPane().add(lbl_r3phase, "cell 2 8,alignx right");

        lbl_r3phaseCalc = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseCalc, "cell 4 8,alignx center");

        lbl_r3phaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3phaseUnit, "cell 5 8,alignx left");

        lbl_r3phaseAct = new JLabel("");
        frame.getContentPane().add(lbl_r3phaseAct, "cell 6 8,alignx right");

        lbl_r3AphaseUnit = new JLabel("radians");
        frame.getContentPane().add(lbl_r3AphaseUnit, "cell 7 8,alignx right");

        btn_plotZ = new JButton("Plot Z Phase Profiles");
        frame.getContentPane().add(btn_plotZ, "cell 8 8");

        btn_estBkgDens = new JButton("<html>Estimate Bkg & &rho;<sub>0</sub></html>");
        frame.getContentPane().add(btn_estBkgDens, "cell 1 10,growx");

        ll_rho0 = new LabeledLabel("<html>&rho;<sub>0</sub> =</html>", "", null);
        frame.getContentPane().add(ll_rho0, "flowx,cell 2 10");

        ll_estBkgPhase = new LabeledLabel("Estimated Background Phase =", "0.0", "radians");
        frame.getContentPane().add(ll_estBkgPhase, "cell 4 10");

        lbl_stepfive = new JLabel("5.");
        frame.getContentPane().add(lbl_stepfive, "cell 0 9");

        btn_calcMagMom = new JButton("Calculate Magnetic Moment");
        frame.getContentPane().add(btn_calcMagMom, "cell 1 9,growx");

        ltf_magMom = new LabeledLabel("|p|=", "0.0", "<html>radians*pixel<sup>3</sup></html>");
        frame.getContentPane().add(ltf_magMom, "cell 2 9");

        btn_removeBkgS5 = new JButton("Remove Bkg");
        frame.getContentPane().add(btn_removeBkgS5, "cell 1 11,growx");

        btn_loadSimImg = new JButton("Load Simulated Images");
        frame.getContentPane().add(btn_loadSimImg, "cell 1 12,growx");

        ltf_snr = new LabeledTextField("SNR=", "0.0", null, 3);
        frame.getContentPane().add(ltf_snr, "cell 1 13");

        ltf_eps12 = new LabeledTextField("<html>&eta;12 =</html>", "0.0", null, 3);
        frame.getContentPane().add(ltf_eps12, "cell 1 14");

        ltf_eps23 = new LabeledTextField("<html>&eta;23 =</html>", "0.0", null, 3);
        frame.getContentPane().add(ltf_eps23, "cell 1 15");

        ll_ReRi = new LabeledLabel("<html>Real(S<sub>i</sub>) =</html>", "", null);
        frame.getContentPane().add(ll_ReRi, "flowx,cell 2 14");

        ll_ImRi = new LabeledLabel("<html>Imag(S<sub>i</sub>) =</html>", "", null);
        frame.getContentPane().add(ll_ImRi, "flowx,cell 2 15");

        ll_momenterror = new LabeledLabel("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&delta;p/p =</html>", "0.0", "%");
        frame.getContentPane().add(ll_momenterror, "cell 1 14");

        ltf_Ri = new LabeledTextField("<html>R<sub>i</sub> =</html>", "0.0", null, 4);
        frame.getContentPane().add(ltf_Ri, "cell 2 13");

        btn_sumRi = new JButton("Sum");
        frame.getContentPane().add(btn_sumRi, "cell 2 13");

        ltf_B0 = new LabeledTextField("B0  =", "0.0", "T", 3);
        frame.getContentPane().add(ltf_B0, "cell 1 16");

        ltf_TELast = new LabeledTextField("    TE_last =", "0.0", "ms", 3);
        frame.getContentPane().add(ltf_TELast, "cell 1 16,alignx right");

        lbl_stepsix = new JLabel("6.");
        frame.getContentPane().add(lbl_stepsix, "cell 0 17");

        btn_loadTE = new JButton("Load First TE Images");
        frame.getContentPane().add(btn_loadTE, "cell 1 17,growx");

        ltf_TEFirst = new LabeledTextField("TE_first =", "0.0", "ms", 4);
        frame.getContentPane().add(ltf_TEFirst, "cell 2 17");

        ll_dChi = new LabeledLabel("<html>&Delta;&Chi; =</html>", "0.0", "ppm");
        frame.getContentPane().add(ll_dChi, "cell 5 17 4 1,alignx left");

        btn_unk = new JButton("TODO: Name");
        frame.getContentPane().add(btn_unk, "cell 1 18,growx");

        ltf_RChi = new LabeledTextField("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;R<sub>&Delta&Chi</sub> =</html>",
                "0.0", "pixels", 4);
        frame.getContentPane().add(ltf_RChi, "cell 2 18");

        ll_a = new LabeledLabel("a =", "0.0", "pixels");
        frame.getContentPane().add(ll_a, "cell 5 18 4 1,alignx left");

        lbl_stepseven = new JLabel("7.");
        frame.getContentPane().add(lbl_stepseven, "cell 0 19");

        btn_loadspinecho = new JButton("Load Spin Echo Images");
        frame.getContentPane().add(btn_loadspinecho, "cell 1 19,growx");

        ltf_sigSE = new LabeledTextField("<html>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&sigma<sub>SE</sub> =</html>",
                "0.0", null, 4);
        frame.getContentPane().add(ltf_sigSE, "cell 2 19");

        /*
         * // "Second Image Center"
         * frame.getContentPane().add(new JLabel("Second Image Center:"),
         * "cell 1 20,alignx right");
         * // "Second Image Center" coordinates
         * ltf_secondImgX = new LabeledTextField("(", "0", null, 4);
         * ltf_secondImgY = new LabeledTextField(",", "0", null, 4);
         * ltf_secondImgZ = new LabeledTextField(",", "0", "-1)", 4);
         * JPanel SECPanel = new JPanel(new MigLayout("", "[][][]", "[]"));
         * SECPanel.add(ltf_secondImgX);
         * SECPanel.add(ltf_secondImgY);
         * SECPanel.add(ltf_secondImgZ);
         * frame.getContentPane().add(SECPanel, "cell 2 20 7 1,alignx left");
         */

        // "V1,SE Region"
        frame.getContentPane().add(new JLabel("<html>V<sub>1,SE</sub> Box Coordinates:</html>"),
                "cell 1 21,alignx right");
        // "V1,SE Region" coordinates
        ltf_v1seX1 = new LabeledTextField("(", "0", null, 4);
        ltf_v1seY1 = new LabeledTextField(",", "0", null, 4);
        ltf_v1seZ1 = new LabeledTextField(",", "0", "-1)", 4);
        ltf_v1seX2 = new LabeledTextField("(", "0", null, 4);
        ltf_v1seY2 = new LabeledTextField(",", "0", null, 4);
        ltf_v1seZ2 = new LabeledTextField(",", "0", "-1)", 4);
        JPanel V1Panel = new JPanel(new MigLayout("", "[][][][][][]", "[]"));
        V1Panel.add(ltf_v1seX1);
        V1Panel.add(ltf_v1seY1);
        V1Panel.add(ltf_v1seZ1);
        V1Panel.add(ltf_v1seX2);
        V1Panel.add(ltf_v1seY2);
        V1Panel.add(ltf_v1seZ2);
        frame.getContentPane().add(V1Panel, "cell 2 21 7 1,alignx left");

        // "V2,SE Region" labels
        frame.getContentPane().add(new JLabel("<html>V<sub>2,SE</sub> Box Coordinates:</html>"),
                "cell 1 22,alignx right");
        // "V2,SE Region" coordinates
        ltf_v2seX1 = new LabeledTextField("(", "0", null, 4);
        ltf_v2seY1 = new LabeledTextField(",", "0", null, 4);
        ltf_v2seZ1 = new LabeledTextField(",", "0", "-1)", 4);
        ltf_v2seX2 = new LabeledTextField("(", "0", null, 4);
        ltf_v2seY2 = new LabeledTextField(",", "0", null, 4);
        ltf_v2seZ2 = new LabeledTextField(",", "0", "-1)", 4);
        JPanel V2Panel = new JPanel(new MigLayout("", "[][][][][][]", "[]"));
        V2Panel.add(ltf_v2seX1);
        V2Panel.add(ltf_v2seY1);
        V2Panel.add(ltf_v2seZ1);
        V2Panel.add(ltf_v2seX2);
        V2Panel.add(ltf_v2seY2);
        V2Panel.add(ltf_v2seZ2);
        frame.getContentPane().add(V2Panel, "cell 2 22 7 1,alignx left");

        btn_estRadSpinEcho = new JButton("Estimate Object Radius From Spin Echo");
        frame.getContentPane().add(btn_estRadSpinEcho, "cell 1 23");

        ll_V0 = new LabeledLabel("<html>V<sub>0</sub> =</html>", "0.0", "<html>pixels<sup>3</sup></html>");
        frame.getContentPane().add(ll_V0, "flowx,cell 2 23");

        ll_rho0SE = new LabeledLabel("<html>&rho<sub>0,SE</sub> =</html>", "0.0", null);
        frame.getContentPane().add(ll_rho0SE, "flowx,cell 4 23 2 1");

        ll_dChiSE = new LabeledLabel("<html>&Delta&Chi =</html>", "0.0", "ppm");
        frame.getContentPane().add(ll_dChiSE, "flowx,cell 2 24");

        ll_aSE = new LabeledLabel("a =", "0.0", "pixels");
        frame.getContentPane().add(ll_aSE, "flowx,cell 4 24");

        ltf_rcx = new LabeledTextField("x=", "0.0", null, 4);
        frame.getContentPane().add(ltf_rcx, "cell 2 1");

        ltf_rcy = new LabeledTextField("y=", "0.0", null, 4);
        frame.getContentPane().add(ltf_rcy, "cell 2 1");

        ltf_rcz = new LabeledTextField("z=", "0.0", "-1 ,", 4);
        frame.getContentPane().add(ltf_rcz, "cell 2 1");

        ltf_M = new LabeledTextField("|M%|=", "50", "%", 2);
        frame.getContentPane().add(ltf_M, "cell 3 1 6 1");

        ldd_MriAxis = new LabeledDropDown<ImageItem.Axis>("<html>B<sub>0</sub> Axis</html>",
                new DefaultComboBoxModel<ImageItem.Axis>(
                        new ImageItem.Axis[] { ImageItem.Axis.X, ImageItem.Axis.Y,
                                ImageItem.Axis.Z }));
        frame.getContentPane().add(ldd_MriAxis, "cell 6 1 2 1,growx");

        // ActionListeners
        btn_loadImages.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.loadStep2Images();
            }
        });
        btn_estCR.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.estimateCenterAndRCenter();
            }
        });
        btn_showRCenter.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.showRCenter();
            }
        });
        btn_genSubpix.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.generateSubpixelImages();
            }
        });
        // btn_removeBkg.addActionListener(new ActionListener() {
        // public void actionPerformed(ActionEvent e) {
        // Calculate_Magnetic_Moment_3D.removeBkgPhase();
        // }
        // });
        btn_removeBkgS5.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.removeBkgPhase();
            }
        });
        btn_estSubC.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.estimateSubpixelCenter();
            }
        });
        btn_redraw.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.redrawCenter();
            }
        });
        btn_verifyRadii.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.verifyRadii();
            }
        });
        btn_plotX.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plotX();
            }
        });
        btn_plotY.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plotY();
            }
        });
        btn_plotZ.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.plotZ();
            }
        });
        btn_estBkgDens.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.estimateBkgAndSpinDensity();
            }
        });
        btn_calcMagMom.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.calculateMagneticMoment();
            }
        });
        btn_loadSimImg.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.loadSimulatedImages();
            }
        });
        btn_sumRi.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.sumRi();
            }
        });
        btn_loadTE.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.loadTEFirstImages();
            }
        });
        btn_unk.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.noname();
            }
        });
        btn_loadspinecho.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.loadSpinEchoImages();
            }
        });
        btn_estRadSpinEcho.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Calculate_Magnetic_Moment_3D.estimateRadiusFromSpinEcho();
            }
        });
    }

}

class LabeledLabel extends JPanel {
    private JLabel IDLabel;
    private JLabel valueLabel;
    private JLabel unitLabel;

    public LabeledLabel(String labelText, String valueText, String unitText) {
        setLayout(new MigLayout("insets 0", "[grow,fill]", "[]"));
        // Create and add the IDLabel component
        IDLabel = new JLabel(labelText);
        add(IDLabel);

        // Create and add the text field component
        valueLabel = new JLabel(valueText);
        add(valueLabel);

        // Create and add the unit IDLabel component (if provided)
        if (unitText != null) {
            unitLabel = new JLabel(unitText);
            if (!valueText.isEmpty()) {
                add(unitLabel);
            }
        }

    }

    // Getter for the text field value
    public String getValue() {
        return valueLabel.getText();
    }

    // Setter for the text field value
    public void setValue(String value) {
        valueLabel.setText(value);
        if (!value.isEmpty() && unitLabel != null)
            add(unitLabel);
    }

    // Additional getters and setters for IDLabel, unit, etc. can be added as needed
}

class LabeledTextField extends JPanel {
    private JLabel label;
    private DefTextField textField;
    private JLabel unitLabel;

    public LabeledTextField(String labelText, String valueText, String unitText, int col) {
        setLayout(new MigLayout("insets 0", "[grow,fill]", "[]"));
        // Create and add the label component
        label = new JLabel(labelText);
        add(label);

        // Create and add the text field component
        textField = (valueText == null) ? (new DefTextField()) : (new DefTextField(valueText));
        textField.setColumns(col);
        add(textField, "growx");

        // Create and add the unit label component (if provided)
        if (unitText != null) {
            unitLabel = new JLabel(unitText);
            add(unitLabel);
        }
    }

    // Getter for the text field
    public DefTextField getTextFieldInstance() {
        return textField;
    }

    // Getter for the text field value
    public String getValue() {
        return textField.getText();
    }

    // Setter for the text field value
    public void setValue(String value) {
        textField.setText(value);
    }

    public void addTFDocumentListener(DocumentListener d) {
        textField.getDocument().addDocumentListener(d);
    }
}

class DefTextField extends JTextField {
    private final String defaultText;

    public DefTextField() {
        super("0.0");
        this.defaultText = "0.0";

        addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                if (getText().isEmpty()) {
                    setText("0.0");
                }
            }
        });
    }

    public DefTextField(String defaultText) {
        super(defaultText);
        this.defaultText = defaultText;

        addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                if (getText().isEmpty()) {
                    setText(defaultText);
                }
            }
        });
    }

    public String getDefaultText() {
        return defaultText;
    }

    public boolean isDefault() {
        return this.getText().compareTo(this.defaultText) == 0;
    }

}

class LabeledDropDown<T> extends JPanel {
    private JLabel label;
    private JComboBox<T> comboBox;

    public LabeledDropDown(String labelText, DefaultComboBoxModel<T> comboBoxModel) {
        setLayout(new MigLayout("insets 0", "[grow,fill]", "[]"));
        // Create and add the label component
        label = new JLabel(labelText);
        add(label);

        // Create and add the combobox component
        comboBox = new JComboBox<T>();
        comboBox.setModel(comboBoxModel);
        add(comboBox);
    }

    // Getter for the text field value
    public T getValue() {
        @SuppressWarnings("unchecked")
        T retval = (T) comboBox.getSelectedItem(); // combobox is uneditable so warning is ok
        return retval;
    }

    // Setter for the text field value
    public void setValue(T value) {
        comboBox.setSelectedItem(value);
    }
}