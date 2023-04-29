import javax.swing.JPanel;
import net.miginfocom.swing.MigLayout;
import javax.swing.JLabel;

public class LabeledLabel extends JPanel {
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
