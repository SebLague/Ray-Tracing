using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class MeshChunk
{
	public Triangle[] triangles;
	public Bounds bounds;
	public int subMeshIndex;

	public MeshChunk(Triangle[] triangles, Bounds bounds, int subMeshIndex)
	{
		this.triangles = triangles;
		this.bounds = bounds;
		this.subMeshIndex = subMeshIndex;
	}
}
