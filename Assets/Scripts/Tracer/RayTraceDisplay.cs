using UnityEngine;

public class RayTraceDisplay : MonoBehaviour
{
	public RayComputeManager raytracer;
	public Shader displayShader;
	Material displayMat;

	void OnRenderImage(RenderTexture source, RenderTexture destination)
	{
		if (raytracer.IsRendering)
		{
			InitMaterial(displayShader, ref displayMat);
			displayMat.SetInt("Frame", raytracer.accumulate ? raytracer.numAccumulatedFrames : 1);

			RenderTexture rt = raytracer.accumulate ? raytracer.accumulatedResult : raytracer.raytraceFrameTex;
			Graphics.Blit(rt, destination, displayMat);
		}
		else
		{
			Graphics.Blit(source, destination);
		}
	}
	
	public static void InitMaterial(Shader shader, ref Material mat)
	{
		if (!mat || (mat.shader != shader && shader))
		{
			if (!shader)
			{
				shader = Shader.Find("Unlit/Texture");
			}

			mat = new Material(shader);
		}
	}
}