module msdf.signeddistance;

import std.algorithm : abs;

class SignedDistance {
public:
    double distance;
    double dot;

    this() {
        distance = double.min_normal;
        dot = 1f;
    }

    this(double dist, double d) {
        distance = dist;
        dot = d;
    }

    int opCmp(const SignedDistance other) const {
        if ((abs(distance) < abs(other.distance)) || (abs(distance) == abs(other.distance) && dot < other.dot))
            return -1;
        else if ((abs(distance) > abs(other.distance)) || (abs(distance) == abs(other.distance) && dot > other.dot))
            return 1;
        else
            return 0;
    }
}