using UnityEngine;

public readonly struct Triangle
{
    public readonly Vector3 PosA;
    public readonly Vector3 PosB;
    public readonly Vector3 PosC;

    public readonly Vector3 NormalA;
    public readonly Vector3 NormalB;
    public readonly Vector3 NormalC;

    public Triangle(Vector3 posA, Vector3 posB, Vector3 posC, Vector3 normalA, Vector3 normalB, Vector3 normalC)
    {
        this.PosA = posA;
        this.PosB = posB;
        this.PosC = posC;
        this.NormalA = normalA;
        this.NormalB = normalB;
        this.NormalC = normalC;
    }
}