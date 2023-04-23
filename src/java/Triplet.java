class Triplet<T> {
    private T first;
    private T second;
    private T third;

    public Triplet(T first, T second, T third){
        this.first = first;
        this.second = second;
        this.third = third;
    }

    public T get(int i){
        switch(i){
            case 0:
                return first;
            case 1:
                return second;
            case 2:
                return third;
            default:
                return null;
        }
    }

    public void set(int i, T val){
        switch (i){
            case 0:
                first = val;
                return;
            case 1:
                second = val;
                return;
            case 2:
                third = val;
                return;
            default:
                return;
        }
    }
}