using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public static class MeshSplitter
{
	const int maxDepth = 6;
	const int maxTrisPerChunk = 48;

	public static MeshChunk[] CreateChunks(Mesh mesh)
	{
		MeshChunk[] subMeshes = new MeshChunk[mesh.subMeshCount];

		Vector3[] verts = mesh.vertices;
		Vector3[] normals = mesh.normals;
		int[] indices = mesh.triangles;

		for (int i = 0; i < subMeshes.Length; i++)
		{
			UnityEngine.Rendering.SubMeshDescriptor subMeshInfo = mesh.GetSubMesh(i);
			var subMeshIndices = indices.AsSpan(subMeshInfo.indexStart, subMeshInfo.indexCount);
			subMeshes[i] = CreateSubMesh(verts, normals, subMeshIndices, i);
		}

		List<MeshChunk> splitChunksList = new List<MeshChunk>();
		foreach (MeshChunk subMesh in subMeshes)
		{
			Split(subMesh, splitChunksList);
		}

		return splitChunksList.ToArray();
	}

	static MeshChunk CreateSubMesh(Vector3[] verts, Vector3[] normals, Span<int> indices, int subMeshIndex)
	{

		Triangle[] triangles = new Triangle[indices.Length / 3];
		Bounds bounds = new Bounds(verts[indices[0]], Vector3.one * 0.01f);

		for (int i = 0; i < indices.Length; i += 3)
		{
			int a = indices[i];
			int b = indices[i + 1];
			int c = indices[i + 2];

			Vector3 posA = verts[a];
			Vector3 posB = verts[b];
			Vector3 posC = verts[c];
			bounds.Encapsulate(posA);
			bounds.Encapsulate(posB);
			bounds.Encapsulate(posC);

			Vector3 normalA = normals[a];
			Vector3 normalB = normals[b];
			Vector3 normalC = normals[c];

			Triangle triangle = new Triangle(posA, posB, posC, normalA, normalB, normalC);
			triangles[i / 3] = triangle;
		}

		return new MeshChunk(triangles, bounds, subMeshIndex);
	}

	static void Split(MeshChunk currChunk, List<MeshChunk> splitChunks, int depth = 0)
	{
		if (currChunk.triangles.Length > maxTrisPerChunk && depth < maxDepth)
		{
			Vector3 q = currChunk.bounds.size / 4;
			Triangle[] allTriangles = currChunk.triangles;
			HashSet<int> takenTriangles = new();

			for (int x = -1; x <= 1; x += 2)
			{
				for (int y = -1; y <= 1; y += 2)
				{
					for (int z = -1; z <= 1; z += 2)
					{
						int remainingTris = allTriangles.Length - takenTriangles.Count;
						if (remainingTris > 0)
						{
							Vector3 splitBoundsOffset = new Vector3(q.x * x, q.y * y, q.z * z);
							Bounds splitBounds = new Bounds(currChunk.bounds.center + splitBoundsOffset, q * 2);

							MeshChunk splitChunk = Extract(allTriangles, takenTriangles, splitBounds, currChunk.subMeshIndex);
							if (splitChunk.triangles.Length > 0)
							{
								Split(splitChunk, splitChunks, depth + 1);
							}
						}
					}
				}
			}
		}
		else
		{
			splitChunks.Add(currChunk);
		}
	}

	static MeshChunk Extract(Triangle[] triangles, HashSet<int> takenTriangles, Bounds splitBounds, int subMeshIndex)
	{
		List<Triangle> newTris = new List<Triangle>();
		Bounds newBounds = new Bounds(splitBounds.center, splitBounds.size);

		for (int i = 0; i < triangles.Length; i++)
		{
			if (takenTriangles.Contains(i))
			{
				continue;
			}

			if (splitBounds.Contains(triangles[i].posA) || splitBounds.Contains(triangles[i].posB) || splitBounds.Contains(triangles[i].posC))
			{
				newBounds.Encapsulate(triangles[i].posA);
				newBounds.Encapsulate(triangles[i].posB);
				newBounds.Encapsulate(triangles[i].posC);
				newTris.Add(triangles[i]);
				takenTriangles.Add(i);
			}
		}

		return new MeshChunk(newTris.ToArray(), newBounds, subMeshIndex);
	}
}
