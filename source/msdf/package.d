/**
    Partial port of msdfgen by Viktor Chlumsky
    by LunaTheFoxgirl and ZyeByte
    Copyright (c) 2023, Inochi2D Project
    
    --- License ---

    MIT License
    Copyright (c) 2016 - 2023 Viktor Chlumsky

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
module msdf;

import msdf.shape;
import msdf.bitmap;
import msdf.contourcombiner;
import msdf.edgeselectors;
import msdf.projection;
import inmath.linalg;

/// Mode of operation.
enum Mode {

    /// Skips error correction pass.
    DISABLED,

    /// Corrects all discontinuities of the distance field regardless if edges are adversely affected.
    INDISCRIMINATE,

    /// Corrects artifacts at edges and other discontinuous distances only if it does not affect edges or corners.
    EDGE_PRIORITY,

    /// Only corrects artifacts at edges.
    EDGE_ONLY
}

/// Configuration of whether to use an algorithm that computes the exact shape distance at the positions of suspected artifacts. This algorithm can be much slower.
enum DistanceCheckMode {

    /// Never computes exact shape distance.
    DO_NOT_CHECK_DISTANCE,

    /// Only computes exact shape distance at edges. Provides a good balance between speed and precision.
    CHECK_DISTANCE_AT_EDGE,

    /// Computes and compares the exact shape distance for each suspected artifact.
    ALWAYS_CHECK_DISTANCE
}

/// The configuration of the MSDF error correction pass.
struct ErrorCorrectionConfig {
    
    /// The default value of minDeviationRatio.
    static  double defaultMinDeviationRatio = 1.11111111111111111;

    /// The default value of minImproveRatio.
    static  double defaultMinImproveRatio = 1.11111111111111111;

    /// Mode of operation.
    Mode mode;
    
    /// Configuration of whether to use an algorithm that computes the exact shape distance at the positions of suspected artifacts. This algorithm can be much slower.
    DistanceCheckMode distanceCheckMode;

    /// The minimum ratio between the actual and maximum expected distance delta to be considered an error.
    double minDeviationRatio;

    /// The minimum ratio between the pre-correction distance error and the post-correction distance error. Has no effect for DO_NOT_CHECK_DISTANCE.
    double minImproveRatio;

    static ErrorCorrectionConfig createDefault() {
        return ErrorCorrectionConfig(Mode.EDGE_PRIORITY, DistanceCheckMode.CHECK_DISTANCE_AT_EDGE, defaultMinDeviationRatio, defaultMinImproveRatio);
    }
}

/// The configuration of the distance field generator algorithm.
struct GeneratorConfig {

    /// Configuration of the error correction pass.
    ErrorCorrectionConfig errorCorrection;

    /// Specifies whether to use the version of the algorithm that supports overlapping contours with the same winding. May be set to false to improve performance when no such contours are present.
    bool overlapSupport;

    this(bool overlapSupport, ErrorCorrectionConfig errorCorrection = ErrorCorrectionConfig.init) {
        this.overlapSupport = overlapSupport;
        this.errorCorrection = errorCorrection;
    }
}

class DistancePixelConversion(T)
    if (is(T == double)) {
private:
    double invRange;

public:
    alias BitmapType = Bitmap!(float, 1);

    void opCall(float* pixels, ref MultiDistance distance) {
        *pixels = cast(float) (invRange * distance * 0.5);
    }
}

class DistancePixelConversion(T) 
    if (is(T == MultiDistance)) {
private:
    double invRange;

public:
    alias BitmapType = Bitmap!(float, 3);

    void opCall(float* pixels, ref MultiDistance distance) {
        pixels[0] = cast(float) (invRange * distance.r + .5);
        pixels[1] = cast(float) (invRange * distance.g + .5);
        pixels[2] = cast(float) (invRange * distance.b + .5);
    }
}

void generateDistanceField(ContourCombiner)(DistancePixelConversion!(ContourCombiner.DistanceType).BitmapType output, Shape shape, Projection projection, double range) {
    auto distancePixelConversion = new DistancePixelConversion!(ContourCombiner.DistanceType)();
    
    auto distanceFinder = new ShapeDistanceFinder!ContourCombiner(shape);
    bool rightToLeft = false;

    for (int y = 0; y < output.height; ++y) {
        int row = shape.inverseYAxis ? output.height-y-1 : y;
        for (int col = 0; col < output.width; ++col) {
            int x = rightToLeft ? output.width-col-1 : col;
            vec2d p = projection.unproject(vec2d(x + 0.5, y + 0.5));
            ContourCombiner.DistanceType distance = distanceFinder.distance(p);
            distancePixelConversion(output(x, row), distance);
        }
        rightToLeft = !rightToLeft;
    }
}

void generateMSDF(ref Bitmap!(float, 3) output, ref Shape shape, ref Projection projection, double range, GeneratorConfig config = GeneratorConfig(true)) {
    if (config.overlapSupport)
        generateDistanceField!(OverlappingContourCombiner!MultiDistanceSelector)(output, shape, projection, range);
    else
        generateDistanceField!(SimpleContourCombiner!MultiDistanceSelector)(output, shape, projection, range);
}