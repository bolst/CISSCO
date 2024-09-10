public enum Axis {
    X, Y, Z, UNKNOWN;

    @Override
    public String toString() {
        switch (this) {
            case X:
                return "X";
            case Y:
                return "Y";
            case Z:
                return "Z";
            case UNKNOWN:
            default:
                return "?";
        }
    }
}