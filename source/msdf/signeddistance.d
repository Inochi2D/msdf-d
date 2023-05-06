module msdf.signeddistance;

import std.math : abs;

struct SignedDistance {
public:
    double distance = double.min_normal;
    double dot = 1f;

    int opCmp(const SignedDistance other) const {
        if ((abs(distance) < abs(other.distance)) || (abs(distance) == abs(other.distance) && dot < other.dot))
            return -1;
        else if ((abs(distance) > abs(other.distance)) || (abs(distance) == abs(other.distance) && dot > other.dot))
            return 1;
        else
            return 0;
    }
}