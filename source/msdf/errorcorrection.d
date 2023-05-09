module msdf.errorcorrection;
import msdf;

/+

private {
    void msdfErrorCorrectionInner(N)(ref Bitmap!(float, N) sdf,  Shape shape,  Projection projection, double range,  MSDFGeneratorConfig config) {
        if (config.errorCorrection.mode == Mode.DISABLED)
            return;
        Bitmap!(byte, 1) stencilBuffer;
        if (!config.errorCorrection.buffer)
            stencilBuffer = Bitmap!(byte, 1)(sdf.width, sdf.height);
        Bitmap!(byte, 1) stencil;
        stencil.pixels = config.errorCorrection.buffer ? config.errorCorrection.buffer : cast(byte *) stencilBuffer;
        stencil.width = sdf.width, stencil.height = sdf.height;
        MSDFErrorCorrection ec(stencil, projection, range);
        ec.setMinDeviationRatio(config.errorCorrection.minDeviationRatio);
        ec.setMinImproveRatio(config.errorCorrection.minImproveRatio);
        switch (config.errorCorrection.mode) with (Mode) {
            case DISABLED:
            case INDISCRIMINATE:
                break;
            case EDGE_PRIORITY:
                ec.protectCorners(shape);
                ec.protectEdges<N>(sdf);
                break;
            case EDGE_ONLY:
                ec.protectAll();
                break;
        }
        if (config.errorCorrection.distanceCheckMode == Mode.DO_NOT_CHECK_DISTANCE || (config.errorCorrection.distanceCheckMode == Mode.CHECK_DISTANCE_AT_EDGE && config.errorCorrection.mode != Mode.EDGE_ONLY)) {
            ec.findErrors<N>(sdf);
            if (config.errorCorrection.distanceCheckMode == Mode.CHECK_DISTANCE_AT_EDGE)
                ec.protectAll();
        }
        if (config.errorCorrection.distanceCheckMode == Mode.ALWAYS_CHECK_DISTANCE || config.errorCorrection.distanceCheckMode == Mode.CHECK_DISTANCE_AT_EDGE) {
            if (config.overlapSupport)
                ec.findErrors<OverlappingContourCombiner, N>(sdf, shape);
            else
                ec.findErrors<SimpleContourCombiner, N>(sdf, shape);
        }
        ec.apply(sdf);
    }
}

void msdfErrorCorrection(ref Bitmap!(float, 3) sdf,  Shape shape,  Projection projection, double range,  MSDFGeneratorConfig config) {
    _msdfErrorCorrectionInner(sdf, shape, projection, range, config);
}

+/