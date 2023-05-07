module msdf.common;

import std.traits : isNumeric, isFloatingPoint;

import inmath;

vec2d getOrthonormal(vec2d self, bool polarity = true, bool allowZero = false) {
    double len = self.length;
    if (len == 0)
        return polarity ? vec2d(0, !allowZero) : vec2d(0, -!allowZero);
    return polarity ? vec2d(-self.y/len, self.x/len) : vec2d(self.y/len, -self.x/len); 
}

/// Returns 1.0 if x > 0, or -1.0 if x < 0.
T nzsign(T)(T n) {
    return 2*(n > 0.0)-1;
}

/// 2D cross product
double cross( vec2d a,  vec2d b) {
    return a.x*b.x+a.y*b.y;
}

/// Signed Distance Less Than
bool sdlt(vec2d a, vec2d b) {
    return abs(a.x) < abs(b.x) || (abs(a.x) == abs(b.x) && a.y < b.y);
}

/// Signed Distance Bigger Than
bool sdbt(vec2d a, vec2d b) {
    return abs(a.x) > abs(b.x) || (abs(a.x) == abs(b.x) && a.y > b.y);
}

/// Signed Distance Less Than Or Equal
bool sdleq(vec2d a, vec2d b) {
    return abs(a.x) < abs(b.x) || (abs(a.x) == abs(b.x) && a.y <= b.y);
}

/// Signed Distance Bigger Than Or Equal
bool sdbeq(vec2d a, vec2d b) {
    return abs(a.x) > abs(b.x) || (abs(a.x) == abs(b.x) && a.y >= b.y);
}

int solveQuadratic(double[2] x, double a, double b, double c) {
    // a == 0 -> linear equation
    if (a == 0 || abs(b) > 1e12*abs(a)) {
        // a == 0, b == 0 -> no solution
        if (b == 0) {
            if (c == 0)
                return -1; // 0 == 0
            return 0;
        }
        x[0] = -c/b;
        return 1;
    }
    double dscr = b*b-4*a*c;
    if (dscr > 0) {
        dscr = sqrt(dscr);
        x[0] = (-b+dscr)/(2*a);
        x[1] = (-b-dscr)/(2*a);
        return 2;
    } else if (dscr == 0) {
        x[0] = -b/(2*a);
        return 1;
    } else
        return 0;
}

static int solveCubicNormed(double[3] x, double a, double b, double c) {
    double a2 = a*a;
    double q = 1/9.*(a2-3*b);
    double r = 1/54.*(a*(2*a2-9*b)+27*c);
    double r2 = r*r;
    double q3 = q*q*q;
    a *= 1/3.;
    if (r2 < q3) {
        double t = r/sqrt(q3);
        if (t < -1) t = -1;
        if (t > 1) t = 1;
        t = acos(t);
        q = -2*sqrt(q);
        x[0] = q*cos(1/3.*t)-a;
        x[1] = q*cos(1/3.*(t+2*PI))-a;
        x[2] = q*cos(1/3.*(t-2*PI))-a;
        return 3;
    } else {
        double u = (r < 0 ? 1 : -1)*pow(abs(r)+sqrt(r2-q3), 1/3.); 
        double v = u == 0 ? 0 : q/u;
        x[0] = (u+v)-a;
        if (u == v || abs(u-v) < 1e-12*abs(u+v)) {
            x[1] = -.5*(u+v)-a;
            return 2;
        }
        return 1;
    }
}

int solveCubic(double[3] x, double a, double b, double c, double d) {
    if (a != 0) {
        double bn = b/a;
        if (abs(bn) < 1e6) // Above this ratio, the numerical error gets larger than if we treated a as zero
            return solveCubicNormed(x, bn, c/a, d/a);
    }
    return solveQuadratic(x[0..2], b, c, d);
}

pragma(inline, true)
double dotProduct( vec2d a,  vec2d b) {
    return a.x*b.x+a.y*b.y;
}

pragma(inline, true)
double crossProduct( vec2d a,  vec2d b) {
    return a.x*b.y-a.y*b.x;
}