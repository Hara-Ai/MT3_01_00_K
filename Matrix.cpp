#include "Matrix.h"
#include <cmath>

Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
    Matrix4x4 matrix = {};
    float tanHalfFovY = tanf(fovY / 2.0f);

    matrix.m[0][0] = 1.0f / (aspectRatio * tanHalfFovY);
    matrix.m[1][1] = 1.0f / tanHalfFovY;
    matrix.m[2][2] = farClip / (farClip - nearClip);
    matrix.m[2][3] = 1.0f;
    matrix.m[3][2] = (-nearClip * farClip) / (farClip - nearClip);

    return matrix;
}

Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
{
    Matrix4x4 matrix = {};

    matrix.m[0][0] = 2.0f / (right - left);
    matrix.m[1][1] = 2.0f / (top - bottom);
    matrix.m[2][2] = 1.0f / (farClip - nearClip);
    matrix.m[3][0] = -(right + left) / (right - left);
    matrix.m[3][1] = -(top + bottom) / (top - bottom);
    matrix.m[3][2] = -nearClip / (farClip - nearClip);
    matrix.m[3][3] = 1.0f;

    return matrix;
}

Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) 
{
    Matrix4x4 matrix = {};

    matrix.m[0][0] = width / 2.0f;
    matrix.m[1][1] = -height / 2.0f;
    matrix.m[2][2] = maxDepth - minDepth;
    matrix.m[3][0] = left + width / 2.0f;
    matrix.m[3][1] = top + height / 2.0f;
    matrix.m[3][2] = minDepth;
    matrix.m[3][3] = 1.0f;

    return matrix;
}


