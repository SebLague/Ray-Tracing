using UnityEngine;
using UnityEngine.Serialization;

[System.Serializable]
public struct RayTracingMaterial
{
	public enum MaterialFlag
	{
		Default,
		CheckerPattern,
		Glass
	}

	[FormerlySerializedAs("colour")]
	public Color diffuseCol;
	[FormerlySerializedAs("emissionColour")]
	public Color emissionCol;
	[FormerlySerializedAs("specularColour")]
	public Color specularCol;

	public Color absorption;
	public float absorptionMultiplier;
	public float emissionStrength;
	[Range(0, 1)] public float smoothness;
	[Range(0, 1)] public float specularProbability;
	public float ior;
	public MaterialFlag flag;

	public void SetDefaultValues()
	{
		diffuseCol = Color.white;
		emissionCol = Color.white;
		emissionStrength = 0;
		specularCol = Color.white;
		smoothness = 0;
		specularProbability = 1;
		ior = 1;
	}
}