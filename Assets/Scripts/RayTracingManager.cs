using UnityEngine;
using UnityEngine.Rendering;
using System.Collections.Generic;

public class RayTracingManager : MonoBehaviour
{
    public enum VisMode
    {
        Default = 0,
        TriangleTestCount = 1,
        BoxTestCount = 2,
        Distance = 3,
        Normal = 4
    }

    [Header("Main Settings")]
    [SerializeField] bool rayTracingEnabled = true;
    public bool accumulate = true;
    public bool useSky;
    [SerializeField] float sunFocus = 500;
    [SerializeField] float sunIntensity = 10;
    [SerializeField] Color sunColor = Color.white;

    [SerializeField, Range(0, 32)] int maxBounceCount = 4;
    [SerializeField, Range(0, 64)] int numRaysPerPixel = 2;
    [SerializeField, Min(0)] float defocusStrength = 0;
    [SerializeField, Min(0)] float divergeStrength = 0.3f;
    [SerializeField, Min(0)] float focusDistance = 1;

    [Header("Debug Settings")]
    [SerializeField] VisMode visMode;
    [SerializeField] float triTestVisScale;
    [SerializeField] float boxTestVisScale;
    [SerializeField] float distanceTestVisScale;
    [SerializeField] bool useSceneView;

    [Header("References")]
    [SerializeField] Shader rayTracingShader;
    [SerializeField] Shader accumulateShader;

    [Header("Info")]
    [SerializeField] int numAccumulatedFrames;

    // Materials and render textures
    Material rayTracingMaterial;
    Material accumulateMaterial;
    RenderTexture resultTexture;

    // Buffers
    ComputeBuffer triangleBuffer;
    ComputeBuffer nodeBuffer;
    ComputeBuffer modelBuffer;

    MeshInfo[] meshInfo;
    Model[] models;
    bool hasBVH;
    LocalKeyword debugVisShaderKeyword;


    private void OnEnable()
    {
        numAccumulatedFrames = 0;
        hasBVH = false;
    }

    private void Update()
    {
        if (Input.GetKeyDown(KeyCode.Space))
        {
            numAccumulatedFrames = 1;
            Debug.Log("Reset render");
        }

        if (Input.GetKeyDown(KeyCode.S))
        {
            string path = System.IO.Path.Combine(Application.persistentDataPath, "screencap_ray.png");
            ScreenCapture.CaptureScreenshot(path);
            Debug.Log("Screenshot: " + path);
        }
    }

    // Called after any camera (e.g. game or scene camera) has finished rendering into the src texture
    void OnRenderImage(RenderTexture src, RenderTexture target)
    {
        if (!Application.isPlaying)
        {
            Graphics.Blit(src, target); // Draw the unaltered camera render to the screen
            return;
        }

        bool isSceneCam = Camera.current.name == "SceneCamera";
        // Debug.Log("Rendering... isscenecam = " + isSceneCam + "  " + Camera.current.name);
        if (isSceneCam)
        {
            if (rayTracingEnabled && useSceneView)
            {
                InitFrame();
                Graphics.Blit(null, target, rayTracingMaterial);
            }
            else
            {
                Graphics.Blit(src, target); // Draw the unaltered camera render to the screen
            }
        }
        else
        {
            Camera.current.cullingMask = rayTracingEnabled ? 0 : 2147483647;
            if (rayTracingEnabled && !useSceneView)
            {
                InitFrame();

                if (accumulate && visMode == VisMode.Default)
                {
                    // Create copy of prev frame
                    RenderTexture prevFrameCopy = RenderTexture.GetTemporary(src.width, src.height, 0, ShaderHelper.RGBA_SFloat);
                    Graphics.Blit(resultTexture, prevFrameCopy);

                    // Run the ray tracing shader and draw the result to a temp texture
                    rayTracingMaterial.SetInt("Frame", numAccumulatedFrames);
                    RenderTexture currentFrame = RenderTexture.GetTemporary(src.width, src.height, 0, ShaderHelper.RGBA_SFloat);
                    Graphics.Blit(null, currentFrame, rayTracingMaterial);

                    // Accumulate
                    accumulateMaterial.SetInt("_Frame", numAccumulatedFrames);
                    accumulateMaterial.SetTexture("_PrevFrame", prevFrameCopy);
                    Graphics.Blit(currentFrame, resultTexture, accumulateMaterial);

                    // Draw result to screen
                    Graphics.Blit(resultTexture, target);

                    // Release temps
                    RenderTexture.ReleaseTemporary(prevFrameCopy);
                    RenderTexture.ReleaseTemporary(currentFrame);
                    numAccumulatedFrames += Application.isPlaying ? 1 : 0;
                }
                else
                {
                    numAccumulatedFrames = 0;
                    Graphics.Blit(null, target, rayTracingMaterial);
                }
            }
            else
            {
                Graphics.Blit(src, target); // Draw the unaltered camera render to the screen
            }
        }
    }

    void InitFrame()
    {
        // Create materials used in blits
        if (rayTracingMaterial == null || rayTracingMaterial.shader != rayTracingShader)
        {
            ShaderHelper.InitMaterial(rayTracingShader, ref rayTracingMaterial);
            debugVisShaderKeyword = new LocalKeyword(rayTracingShader, "DEBUG_VIS");
        }
        ShaderHelper.InitMaterial(accumulateShader, ref accumulateMaterial);
        ShaderHelper.CreateRenderTexture(ref resultTexture, Screen.width, Screen.height, FilterMode.Bilinear, ShaderHelper.RGBA_SFloat, "Result");
        models = FindObjectsOfType<Model>();

        if (!hasBVH)
        {
            var data = CreateAllMeshData(models);
            hasBVH = true;

            meshInfo = data.meshInfo.ToArray();
            ShaderHelper.CreateStructuredBuffer(ref modelBuffer, meshInfo);

            // Triangles buffer
            ShaderHelper.CreateStructuredBuffer(ref triangleBuffer, data.triangles);
            rayTracingMaterial.SetBuffer("Triangles", triangleBuffer);
            rayTracingMaterial.SetInt("triangleCount", triangleBuffer.count);

            // Node buffer
            ShaderHelper.CreateStructuredBuffer(ref nodeBuffer, data.nodes);
            rayTracingMaterial.SetBuffer("Nodes", nodeBuffer);
        }
        UpdateModels();
        // Update data
        UpdateCameraParams(Camera.current);
        SetShaderParams();
    }

    void SetShaderParams()
    {
        rayTracingMaterial.SetKeyword(debugVisShaderKeyword, visMode != VisMode.Default);
        rayTracingMaterial.SetInt("visMode", (int)visMode);
        float debugVisScale = visMode switch
        {
            VisMode.TriangleTestCount => triTestVisScale,
            VisMode.BoxTestCount => boxTestVisScale,
            VisMode.Distance => distanceTestVisScale,
            _ => triTestVisScale
        };
        rayTracingMaterial.SetFloat("debugVisScale", debugVisScale);
        rayTracingMaterial.SetInt("Frame", numAccumulatedFrames);
        rayTracingMaterial.SetInt("UseSky", useSky ? 1 : 0);

        rayTracingMaterial.SetInt("MaxBounceCount", maxBounceCount);
        rayTracingMaterial.SetInt("NumRaysPerPixel", numRaysPerPixel);
        rayTracingMaterial.SetFloat("DefocusStrength", defocusStrength);
        rayTracingMaterial.SetFloat("DivergeStrength", divergeStrength);

        rayTracingMaterial.SetFloat("SunFocus", sunFocus);
        rayTracingMaterial.SetFloat("SunIntensity", sunIntensity);
        rayTracingMaterial.SetColor("SunColour", sunColor);
    }

    void UpdateCameraParams(Camera cam)
    {
        float planeHeight = focusDistance * Mathf.Tan(cam.fieldOfView * 0.5f * Mathf.Deg2Rad) * 2;
        float planeWidth = planeHeight * cam.aspect;
        // Send data to shader
        rayTracingMaterial.SetVector("ViewParams", new Vector3(planeWidth, planeHeight, focusDistance));
        rayTracingMaterial.SetMatrix("CamLocalToWorldMatrix", cam.transform.localToWorldMatrix);
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
        rayTracingMaterial.SetBuffer("ModelInfo", modelBuffer);
        rayTracingMaterial.SetInt("modelCount", models.Length);
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

                BVH bvh = new(model.Mesh.vertices, model.Mesh.triangles, model.Mesh.normals);
                if (model.logBVHStats) Debug.Log($"BVH Stats: {model.gameObject.name}\n{bvh.stats}");

                allData.triangles.AddRange(bvh.GetTriangles());
                allData.nodes.AddRange(bvh.GetNodes());
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

    class MeshDataLists
    {

        public List<Triangle> triangles = new();
        public List<BVH.Node> nodes = new();
        public List<MeshInfo> meshInfo = new();
    }

    void OnDestroy()
    {
        if (Application.isPlaying)
        {
            ShaderHelper.Release(triangleBuffer, nodeBuffer, modelBuffer);
            ShaderHelper.Release(resultTexture);
            Destroy(rayTracingMaterial);
        }
    }

    void OnValidate()
    {
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
