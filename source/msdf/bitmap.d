module msdf.bitmap;

class Bitmap(T, int N = 1) {
private:
    T[] pixels;
    int w, h;

public:
    this() {
        pixels = null;
        w = 0;
        h = 0;
    }

    this(int width, int height) {
        w = width;
        h = height;

        pixels = new T[N*w*h];
    }

    this(Bitmap!(T, N) orig) {
        w = orig.w;
        h = orig.h;

        pixels = orig.pixels.dup;
    }

    version(none)
    Bitmap!(T, N) opAssign(Bitmap!(T, N) orig) {
        if (pixels != orig.pixels) {
            w = orig.w;
            h = orig.h;

            pixels = orig.pixels.dup;
        }

        return this;
    }

    int width()  {
        return w;
    }

    int height()  {
        return h;
    }

    T* opCall(int x, int y) {
        return &pixels[N * (w*y+x)];
    }

    T* opCast(T)() {
        return pixels.ptr;
    }
}