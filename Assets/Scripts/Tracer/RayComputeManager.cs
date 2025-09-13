using UnityEngine;
using System.Collections.Generic;
using Seb.AccelerationStructures;
using Seb.Helpers;
using UnityEngine.Experimental.Rendering;

public class RayComputeManager : MonoBehaviour
{
	[Header("Main Settings")]
	public bool rayTracingEnabled = true;

	public bool accumulate = true;
	public BVH.Quality bvhQuality = BVH.Quality.High;

	[SerializeField, Range(0, 32)] int maxBounceCount = 4;
	[SerializeField] int numRaysPerPixel = 1;
	[SerializeField, Min(0)] float defocusStrength = 0;
	[SerializeField, Min(0)] float divergeStrength = 0.3f;
	[Min(0)] public float focusDistance = 1;

	[Header("Sky Settings")]
	public bool useSky;

	[SerializeField] float sunFocus = 500;
	[SerializeField] float sunIntensity = 10;
	[SerializeField] Color sunColor = Color.white;
	public Transform sunTransform;

	[Header("Debug Settings")]
	public string screenshotName = "rayScreenshot";
	public Vector4 debugParams;

	[Header("References")]
	[SerializeField] ComputeShader rayComputeShader;

	[Header("Info")]
	public int numAccumulatedFrames;

	public int renderSeed;

	[SerializeField] float info_timeSinceReset;
	[SerializeField] Vector2Int screenSize;

	// Buffers
	ComputeBuffer triangleBuffer;
	ComputeBuffer nodeBuffer;
	ComputeBuffer modelBuffer;

	MeshInfo[] meshInfo;
	Model[] models;
	bool hasBVH;

	[HideInInspector] public RenderTexture raytraceFrameTex;
	[HideInInspector] public RenderTexture accumulatedResult;
	public bool IsRendering => Application.isPlaying && rayTracingEnabled;

	const int kernelRayTrace = 0;
	const int kernelResetAccumulated = 1;


	private void OnEnable()
	{
		hasBVH = false;
		renderSeed = new System.Random().Next();

		ResetAccumulatedRender();
	}

	public void ResetAccumulatedRender()
	{
		numAccumulatedFrames = 1;
		info_timeSinceReset = 0;
		InitFrame();

		Seb.Helpers.ComputeHelper.Dispatch(rayComputeShader, accumulatedResult.width, accumulatedResult.height, kernelIndex: kernelResetAccumulated);
	}

	private void Update()
	{
		RenderFrame();
		HandleInput();
	}

	void RenderFrame()
	{
		if (!IsRendering) return;

		InitFrame();

		Seb.Helpers.ComputeHelper.Dispatch(rayComputeShader, raytraceFrameTex.width, raytraceFrameTex.height, kernelIndex: 0);

		info_timeSinceReset += Time.deltaTime;

		if (accumulate) numAccumulatedFrames++;
	}


	void HandleInput()
	{
		if (Input.GetKeyDown(KeyCode.Space))
		{
			ResetAccumulatedRender();
			Debug.Log("Reset render");
		}

		if (Input.GetKeyDown(KeyCode.S))
		{
			string path = System.IO.Path.Combine(Application.persistentDataPath, screenshotName + ".png");
			ScreenCapture.CaptureScreenshot(path);
			Debug.Log("Screenshot: " + path);
		}
	}
	

	void InitFrame()
	{
		InitTexturesAndBuffers();
		models = FindObjectsByType<Model>(FindObjectsInactive.Exclude, FindObjectsSortMode.InstanceID);

		InitBVH();
		UpdateModels();
		UpdateCameraParams(Camera.main);
		SetShaderParams();
	}

	void InitTexturesAndBuffers()
	{
		int width = Screen.width;
		int height = Screen.height;
		screenSize = new Vector2Int(width, height);

		Seb.Helpers.ComputeHelper.CreateRenderTexture(ref raytraceFrameTex, width, height, FilterMode.Bilinear, GraphicsFormat.R32G32B32A32_SFloat, "Raytrace Frame");
		Seb.Helpers.ComputeHelper.CreateRenderTexture(ref accumulatedResult, width, height, FilterMode.Bilinear, GraphicsFormat.R32G32B32A32_SFloat, "Raytrace Accumulated");

		rayComputeShader.SetTexture(kernelRayTrace, "FrameRender", raytraceFrameTex);
		rayComputeShader.SetTexture(kernelRayTrace, "AccumulatedRender", accumulatedResult);
		rayComputeShader.SetTexture(kernelResetAccumulated, "AccumulatedRender", accumulatedResult);

		rayComputeShader.SetInts("Resolution", raytraceFrameTex.width, raytraceFrameTex.height);
		rayComputeShader.SetVector("debugParams", debugParams);
	}

	void InitBVH()
	{
		if (hasBVH) return;

		hasBVH = true;
		var data = CreateAllMeshData(models);

		meshInfo = data.meshInfo.ToArray();
		ComputeHelper.CreateStructuredBuffer(ref modelBuffer, meshInfo);

		// Triangles buffer
		ComputeHelper.CreateStructuredBuffer(ref triangleBuffer, data.triangles);
		rayComputeShader.SetBuffer(kernelRayTrace, "Triangles", triangleBuffer);
		rayComputeShader.SetInt("triangleCount", triangleBuffer.count);

		// Node buffer
		ComputeHelper.CreateStructuredBuffer(ref nodeBuffer, data.nodes);
		rayComputeShader.SetBuffer(kernelRayTrace, "Nodes", nodeBuffer);
	}

	void SetShaderParams()
	{
		rayComputeShader.SetInt("Frame", numAccumulatedFrames);
		rayComputeShader.SetInt("UseSky", useSky ? 1 : 0);

		rayComputeShader.SetInt("MaxBounceCount", maxBounceCount);
		rayComputeShader.SetInt("NumRaysPerPixel", numRaysPerPixel);
		rayComputeShader.SetFloat("DefocusStrength", defocusStrength);
		rayComputeShader.SetFloat("DivergeStrength", divergeStrength);

		rayComputeShader.SetFloat("SunFocus", sunFocus);
		rayComputeShader.SetFloat("SunIntensity", sunIntensity);
		rayComputeShader.SetVector("SunColour", sunColor);
		rayComputeShader.SetVector("dirToSun", sunTransform == null ? Vector3.down : -sunTransform.forward);

		rayComputeShader.SetInt("Frame", numAccumulatedFrames);
		rayComputeShader.SetInt("renderSeed", renderSeed);
		rayComputeShader.SetBool("accumulate", accumulate);
	}

	void UpdateCameraParams(Camera cam)
	{
		float planeHeight = focusDistance * Mathf.Tan(cam.fieldOfView * 0.5f * Mathf.Deg2Rad) * 2;
		float planeWidth = planeHeight * cam.aspect;
		// Send data to shader
		rayComputeShader.SetVector("ViewParams", new Vector3(planeWidth, planeHeight, focusDistance));
		rayComputeShader.SetMatrix("CamLocalToWorldMatrix", cam.transform.localToWorldMatrix);
	}

	void UpdateModels()
	{
		for (int i = 0; i < models.Length; i++)
		{
			meshInfo[i].WorldToLocalMatrix = models[i].transform.worldToLocalMatrix;
			meshInfo[i].LocalToWorldMatrix = models[i].transform.localToWorldMatrix;
			meshInfo[i].Material = models[i].material;
		}

		modelBuffer.SetData(meshInfo);
		rayComputeShader.SetBuffer(kernelRayTrace, "ModelInfo", modelBuffer);
		rayComputeShader.SetInt("modelCount", models.Length);
	}

	MeshDataLists CreateAllMeshData(Model[] models)
	{
		MeshDataLists allData = new();
		Dictionary<Mesh, (int nodeOffset, int triOffset)> meshLookup = new();

		foreach (Model model in models)
		{
			// Construct BVH if this is the first time seeing the current mesh (otherwise reuse)
			if (!meshLookup.ContainsKey(model.Mesh))
			{
				meshLookup.Add(model.Mesh, (allData.nodes.Count, allData.triangles.Count));

				BVH bvh = new(model.Mesh.vertices, model.Mesh.triangles, model.Mesh.normals, bvhQuality);
				if (model.logBVHStats) Debug.Log($"BVH Stats: {model.gameObject.name}\n{bvh.stats}");

				allData.triangles.AddRange(bvh.Triangles);
				allData.nodes.AddRange(bvh.Nodes);
			}

			// Create the mesh info
			allData.meshInfo.Add(new MeshInfo()
			{
				NodeOffset = meshLookup[model.Mesh].nodeOffset,
				TriangleOffset = meshLookup[model.Mesh].triOffset,
				WorldToLocalMatrix = model.transform.worldToLocalMatrix,
				Material = model.material
			});
		}

		return allData;
	}

	void OnDestroy()
	{
		if (Application.isPlaying)
		{
			ComputeHelper.Release(triangleBuffer, nodeBuffer, modelBuffer);
			ComputeHelper.Release(accumulatedResult, raytraceFrameTex);
		}

		Seb.Helpers.ComputeHelper.Release(raytraceFrameTex);
	}

	class MeshDataLists
	{
		public List<BVH.Triangle> triangles = new();
		public List<BVH.Node> nodes = new();
		public List<MeshInfo> meshInfo = new();
	}

	struct MeshInfo
	{
		public int NodeOffset;
		public int TriangleOffset;
		public Matrix4x4 WorldToLocalMatrix;
		public Matrix4x4 LocalToWorldMatrix;
		public RayTracingMaterial Material;
	}
}