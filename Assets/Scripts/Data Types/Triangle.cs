using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public struct Triangle
{
	public Vector3 posA;
	public Vector3 posB;
	public Vector3 posC;

	public Vector3 normalA;
	public Vector3 normalB;
	public Vector3 normalC;

	public Triangle(Vector3 posA, Vector3 posB, Vector3 posC, Vector3 normalA, Vector3 normalB, Vector3 normalC)
	{
		this.posA = posA;
		this.posB = posB;
		this.posC = posC;
		this.normalA = normalA;
		this.normalB = normalB;
		this.normalC = normalC;
	}
}