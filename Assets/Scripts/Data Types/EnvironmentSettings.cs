using UnityEngine;

[System.Serializable]
public struct EnvironmentSettings
{
	public bool enabled;
	public Color groundColour;
	public Color skyColourHorizon;
	public Color skyColourZenith;
	public float sunFocus;
	public float sunIntensity;
}
