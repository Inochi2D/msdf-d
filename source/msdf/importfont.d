module msdf.importfont;

import std.string : toStringz;
import core.stdc.stdint;
import core.stdc.string;

import inmath.linalg;
import bindbc.freetype;

import msdf.shape;
import msdf.contour;
import msdf.segment;

alias unicode_t = uint;

struct GlyphIndex {
private:
    uint index;

public:
    uint getIndex() const {
        return index;
    }
}

struct FontMetrics {
    /// The size of one EM.
    double emSize;
    /// The vertical position of the ascender and descender relative to the baseline.
    double ascenderY, descenderY;
    /// The vertical difference between consecutive baselines.
    double lineHeight;
    /// The vertical position and thickness of the underline.
    double underlineY, underlineThickness;
}

struct FontVariationAxis {
    /// The name of the variation axis.
    const(char)* name;
    /// The axis's minimum coordinate value.
    double minValue;
    /// The axis's maximum coordinate value.
    double maxValue;
    /// The axis's default coordinate value. FreeType computes meaningful default values for Adobe MM fonts.
    double defaultValue;
}

struct FreetypeHandle {
    FT_Library library;
}

struct FontHandle {
    FT_Face face;
    bool ownership;
}

struct FtContext {
    vec2d position;
    Shape shape;
    Contour* contour;
}

private extern(C) nothrow
{
    pragma(inline, true)
    double F26DOT6_TO_DOUBLE(T)(T x) {
        return 1 / 64.0 * cast(double) x;
    }

    vec2d _msdfFtPoint2(FT_Vector vector) {
        return vec2d(F26DOT6_TO_DOUBLE(vector.x), F26DOT6_TO_DOUBLE(vector.y));
    }

    int _msdfFtMoveTo(const(FT_Vector)* to, void* user) {
        FtContext* context = cast(FtContext*) user;
        if (!(context.contour && context.contour.edges.length > 0))
            context.contour = context.shape.addContour();
        context.position = _msdfFtPoint2(*to);
        return 0;
    }

    int _msdfFtLineTo(const(FT_Vector)* to, void* user) {
        FtContext* context = cast(FtContext*) user;
        vec2d endpoint = _msdfFtPoint2(*to);
        if (endpoint != context.position) {
            context.contour.addEdge(new LinearSegment(context.position, endpoint));
            context.position = endpoint;
        }
        return 0;
    }

    int _msdfFtConicTo(const(FT_Vector)* control, const(FT_Vector)* to, void* user) {
        FtContext* context = cast(FtContext*) user;
        context.contour.addEdge(new QuadraticSegment(context.position, _msdfFtPoint2(*control), _msdfFtPoint2(*to)));
        context.position = _msdfFtPoint2(*to);
        return 0;
    }

    int _msdfFtCubicTo(const(FT_Vector)* control1, const(FT_Vector)* control2, const(FT_Vector)* to, void* user) {
        FtContext* context = cast(FtContext*) user;
        context.contour.addEdge(new CubicSegment(context.position, _msdfFtPoint2(*control1), _msdfFtPoint2(*control2), _msdfFtPoint2(*to)));
        context.position = _msdfFtPoint2(*to);
        return 0;
    }
}

FreetypeHandle* initializeFreetype() {
    FreetypeHandle* handle = new FreetypeHandle();
    FT_Error error = FT_Init_FreeType(&handle.library);
    if (error)
        return null;

    return handle;
}

void deinitializeFreetype(FreetypeHandle* library) {
    FT_Done_FreeType(library.library);
}

FontHandle* adoptFreetypeFont(FT_Face ftFace) {
    FontHandle* handle = new FontHandle();
    handle.face = ftFace;
    handle.ownership = false;
    return handle;
}

FT_Error readFreetypeOutline(Shape output, FT_Outline* outline) {
    output.contours.length = 0;
    output.inverseYAxis = false;
    FtContext context;

    context.shape = output;
    FT_Outline_Funcs ftFunctions;

    ftFunctions.move_to = &_msdfFtMoveTo;
    ftFunctions.line_to = &_msdfFtLineTo;
    ftFunctions.conic_to = &_msdfFtConicTo;
    ftFunctions.cubic_to = &_msdfFtCubicTo;
    ftFunctions.shift = 0;
    ftFunctions.delta = 0;

    FT_Error error = FT_Outline_Decompose(outline, &ftFunctions, &context);
    if (output.contours.length > 0 && output.contours[$ - 1].edges.length == 0)
        output.contours.length -= 1;

    return error;
}

FontHandle* loadFont(FreetypeHandle* library, string filename) {
    if (!library)
        return null;

    FontHandle* handle = new FontHandle();
    FT_Error error = FT_New_Face(library.library, filename.toStringz, 0, &handle.face);
    if (error)
        return null;

    handle.ownership = true;
    return handle;
}

FontHandle* loadFontData(FreetypeHandle* library, const(ubyte)[] data) {
    if (!library)
        return null;

    FontHandle* handle = new FontHandle();
    FT_Error error = FT_New_Memory_Face(library.library, data.ptr, cast(int) data.length, 0, &handle.face);
    if (error)
        return null;
    handle.ownership = true;
    return handle;
}

void destroyFont(FontHandle* font) {
    if (font.ownership)
        FT_Done_Face(font.face);
}

bool getFontMetrics(ref FontMetrics metrics, FontHandle* font) {
    metrics.emSize = F26DOT6_TO_DOUBLE(font.face.units_per_EM);
    metrics.ascenderY = F26DOT6_TO_DOUBLE(font.face.ascender);
    metrics.descenderY = F26DOT6_TO_DOUBLE(font.face.descender);
    metrics.lineHeight = F26DOT6_TO_DOUBLE(font.face.height);
    metrics.underlineY = F26DOT6_TO_DOUBLE(font.face.underline_position);
    metrics.underlineThickness = F26DOT6_TO_DOUBLE(font.face.underline_thickness);
    return true;
}

bool getFontWhitespaceWidth(double spaceAdvance, double tabAdvance, FontHandle* font) {
    FT_Error error = FT_Load_Char(font.face, ' ', FT_LOAD_NO_SCALE);
    if (error)
        return false;
    spaceAdvance = F26DOT6_TO_DOUBLE(font.face.glyph.advance.x);
    error = FT_Load_Char(font.face, '\t', FT_LOAD_NO_SCALE);
    if (error)
        return false;
    tabAdvance = F26DOT6_TO_DOUBLE(font.face.glyph.advance.x);
    return true;
}

bool getGlyphIndex(ref GlyphIndex glyphIndex, FontHandle* font, unicode_t unicode) {
    glyphIndex = GlyphIndex(FT_Get_Char_Index(font.face, unicode));
    return glyphIndex.getIndex() != 0;
}

bool loadGlyph(ref Shape output, FontHandle* font, GlyphIndex glyphIndex, double* advance = null) {
    if (!font)
        return false;

    FT_Error error = FT_Load_Glyph(font.face, glyphIndex.getIndex(), FT_LOAD_NO_SCALE);
    if (error)
        return false;

    if (advance)
        *advance = F26DOT6_TO_DOUBLE(font.face.glyph.advance.x);
    
    return !readFreetypeOutline(output, &font.face.glyph.outline);
}

bool loadGlyph(ref Shape output, FontHandle* font, unicode_t unicode, double* advance = null) {
    return loadGlyph(output, font, GlyphIndex(FT_Get_Char_Index(font.face, unicode)), advance);
}

bool getKerning(ref double output, FontHandle* font, GlyphIndex glyphIndex1, GlyphIndex glyphIndex2) {
    // Based on the FT_Kerning_Mode enum, FT_KERNING_UNSCALED is 2.
    enum FT_KERNING_UNSCALED = 2;

    FT_Vector kerning;
    if (FT_Get_Kerning(font.face, glyphIndex1.getIndex(), glyphIndex2.getIndex(), FT_KERNING_UNSCALED, &kerning)) {
        output = 0;
        return false;
    }

    output = F26DOT6_TO_DOUBLE(kerning.x);
    return true;
}

bool getKerning(ref double output, FontHandle* font, unicode_t unicode1, unicode_t unicode2) {
    return getKerning(output, font, GlyphIndex(FT_Get_Char_Index(font.face, unicode1)), GlyphIndex(FT_Get_Char_Index(font.face, unicode2)));
}

bool setFontVariationAxis(FreetypeHandle* library, FontHandle* font, string name, double coordinate) {
    bool success = false;

    if (font.face.face_flags & FT_FACE_FLAG_MULTIPLE_MASTERS) {
        FT_MM_Var* master = null;
        if (FT_Get_MM_Var(font.face, &master))
            return false;

        if (master && master.num_axis > 0) {
            FT_Fixed[] coords;
            coords.length = master.num_axis;

            if (!FT_Get_Var_Design_Coordinates(font.face, cast(FT_UInt) coords.length, &coords[0])) {
                const char* zName = name.toStringz;
                for (FT_UInt i = 0; i < master.num_axis; ++i) {
                    if (!strcmp(zName, master.axis[i].name)) {
                        coords[i] = cast(FT_Fixed) (65536.0 * coordinate);
                        success = true;
                        break;
                    }
                }
            }

            if (FT_Set_Var_Design_Coordinates(font.face, cast(FT_UInt) coords.length, &coords[0]))
                success = false;
        }
        FT_Done_MM_Var(library.library, master);
    }

    return success;
}

bool listFontVariationAxes(FontVariationAxis[] axes, FreetypeHandle* library, FontHandle* font) {
    if (font.face.face_flags & FT_FACE_FLAG_MULTIPLE_MASTERS) {
        FT_MM_Var* master = null;
        if (FT_Get_MM_Var(font.face, &master))
            return false;

        axes.length = master.num_axis;

        for (FT_UInt i = 0; i < master.num_axis; ++i) {
            FontVariationAxis* axis = &axes[i];
            axis.name = master.axis[i].name;
            axis.minValue = master.axis[i].minimum;
            axis.maxValue = master.axis[i].maximum;
            axis.defaultValue = master.axis[i].def;
        }
        FT_Done_MM_Var(library.library, master);
        return true;
    }
    return false;
}