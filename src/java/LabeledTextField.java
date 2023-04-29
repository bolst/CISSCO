import javax.swing.JPanel;
import net.miginfocom.swing.MigLayout;
import javax.swing.JLabel;
import javax.swing.JTextField;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

public class LabeledTextField extends JPanel {
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
    public DefTextField getValueTF() {
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

    // Additional getters and setters for label, unit, etc. can be added as needed
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
