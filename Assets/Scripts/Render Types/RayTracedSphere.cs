using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RayTracedSphere : MonoBehaviour
{
	public RayTracingMaterial material;

	[SerializeField, HideInInspector] int materialObjectID;
	[SerializeField, HideInInspector] bool materialInitFlag;

	void OnValidate()
	{
		if (!materialInitFlag)
		{
			materialInitFlag = true;
			material.SetDefaultValues();
		}

		MeshRenderer renderer = GetComponent<MeshRenderer>();
		if (renderer != null)
		{
			if (materialObjectID != gameObject.GetInstanceID())
			{
				renderer.sharedMaterial = new Material(renderer.sharedMaterial);
				materialObjectID = gameObject.GetInstanceID();
			}
			renderer.sharedMaterial.color = material.colour;
		}
	}
}
