module msdf.projection;

import inmath;

class Projection {
public:
    this() {
        scale = vec2d(1);
        translate = vec2d(0);
    }

    this(in vec2d scale, in vec2d translate) {
        this.scale = scale;
        this.translate = translate;
    }

    vec2d project(in vec2d coord) const {
        return scale * (coord + translate);
    }

    vec2d unproject(in vec2d coord) const {
        return vec2d(coord.x / scale.x, coord.y / scale.y) - translate;
    }

    vec2d projectVector(in vec2d vector) const {
        return scale * vector;
    }

    vec2d unprojectVector(in vec2d vector) const {
        return vec2d(vector.x / scale.x, vector.y - scale.y);
    }

    double projectX(double x) const {
        return scale.x * (x + translate.x);
    }

    double projectY(double y) const {
        return scale.y * (y + translate.y);
    }

    double unprojectX(double x) const {
        return x / scale.x - translate.x;
    }

    double unprojectY(double y) const {
        return y / scale.y - translate.y;
    }

private:
    vec2d scale;
    vec2d translate;
}