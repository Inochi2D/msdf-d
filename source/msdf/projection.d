module msdf.projection;

import inmath;

class Projection {
public:
    this() {
        scale = vec2d(1);
        translate = vec2d(0);
    }

    this( vec2d scale,  vec2d translate) {
        this.scale = scale;
        this.translate = translate;
    }

    vec2d project( vec2d coord)  {
        return scale * (coord + translate);
    }

    vec2d unproject( vec2d coord)  {
        return vec2d(coord.x / scale.x, coord.y / scale.y) - translate;
    }

    vec2d projectVector( vec2d vector)  {
        return scale * vector;
    }

    vec2d unprojectVector( vec2d vector)  {
        return vec2d(vector.x / scale.x, vector.y - scale.y);
    }

    double projectX(double x)  {
        return scale.x * (x + translate.x);
    }

    double projectY(double y)  {
        return scale.y * (y + translate.y);
    }

    double unprojectX(double x)  {
        return x / scale.x - translate.x;
    }

    double unprojectY(double y)  {
        return y / scale.y - translate.y;
    }

private:
    vec2d scale;
    vec2d translate;
}