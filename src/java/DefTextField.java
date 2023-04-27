import javax.swing.JTextField;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

public class DefTextField extends JTextField {
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
