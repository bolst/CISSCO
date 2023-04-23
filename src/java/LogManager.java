import java.io.File;
import java.io.FileWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class LogManager {

    private static File logFile;
    private static ArrayList<String> warnings;
    private static ArrayList<String> info;
    private static ArrayList<String> errors;
    private static ArrayList<String> messages;
    private static String time;
    private static final String WARNINGLABEL = "<WARNING>";
    private static final String INFOLABEL = "<INFO>";
    private static final String ERRORLABEL = "<ERROR>";
    private static final String VARIABLELABEL = "<VARIABLE>";
    private static String currentClass;

    public LogManager() {
        warnings = new ArrayList<String>();
        info = new ArrayList<String>();
        errors = new ArrayList<String>();
        messages = new ArrayList<String>();

        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        LocalDateTime now = LocalDateTime.now();
        time = "(" + dtf.format(now) + ")";
        logFile = new File("log.txt");
        messages.add("--------------------------------------------------------------\n" + "LOG FILE\n" + time
                + "\n--------------------------------------------------------------\n");
        final StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        currentClass = "<" + ste[ste.length - 1].getMethodName() + ">";
    }

    public LogManager(String s) {
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        LocalDateTime now = LocalDateTime.now();
        time = dtf.format(now);
        logFile = new File(s);
    }

    private static void updateCurrentClassAndTime() {
        final StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        currentClass = "<" + ste[ste.length - 1].getMethodName() + ">";
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        LocalDateTime now = LocalDateTime.now();
        time = "(" + dtf.format(now) + ")";
    }

    public void addWarning(String s) {
        updateCurrentClassAndTime();
        warnings.add(time + currentClass + WARNINGLABEL + s);
        messages.add(time + currentClass + WARNINGLABEL + s);
        writeToLogFile();
    }

    public void addInfo(String s) {
        updateCurrentClassAndTime();
        info.add(time + currentClass + INFOLABEL + s);
        messages.add(time + currentClass + INFOLABEL + s);
        writeToLogFile();
    }

    public void addInfo(String s, Object o) {
        updateCurrentClassAndTime();
        info.add(time + currentClass + INFOLABEL + s + ": " + o.toString());
        messages.add(time + currentClass + INFOLABEL + s + ": " + o.toString());
        writeToLogFile();
    }

    public void addVariable(String s, Object o) {
        updateCurrentClassAndTime();
        info.add(time + currentClass + VARIABLELABEL + s + ": " + o.toString());
        messages.add(time + currentClass + VARIABLELABEL + s + ": " + o.toString());
        writeToLogFile();
    }

    public void addError(String s) {
        updateCurrentClassAndTime();
        errors.add(time + currentClass + ERRORLABEL + s);
        messages.add(time + currentClass + ERRORLABEL + s);
        writeToLogFile();
    }

    private static int writeToLogFile() {
        try (FileWriter fw = new FileWriter(logFile)) {
            for (String s : messages) {
                fw.write(s + "\n");
            }
        } catch (Exception exc) {
            return -1;
        }

        return 0;
    }
}
