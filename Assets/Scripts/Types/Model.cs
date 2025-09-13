using UnityEngine;

public class Model : MonoBehaviour
{
	public RayTracingMaterial material;

	[Header("References")]
	public MeshFilter meshFilter;

	public MeshRenderer meshRenderer;
	public bool logBVHStats;
	[SerializeField, HideInInspector] int materialObjectID;

	public Mesh Mesh => meshFilter.sharedMesh;

	private void OnValidate()
	{
		meshFilter = GetComponent<MeshFilter>();
		meshRenderer = GetComponent<MeshRenderer>();
		SetUpMaterialDisplay();
	}

	void SetUpMaterialDisplay()
	{
		if (gameObject.GetInstanceID() != materialObjectID)
		{
			materialObjectID = gameObject.GetInstanceID();
			Material[] originalMaterials = meshRenderer.sharedMaterials;
			Material[] newMaterials = new Material[originalMaterials.Length];
			Shader shader = Shader.Find("Standard");
			for (int i = 0; i < meshRenderer.sharedMaterials.Length; i++)
			{
				newMaterials[i] = new Material(shader);
			}

			meshRenderer.sharedMaterials = newMaterials;
		}


		RayTracingMaterial mat = material;
		bool displayEmissiveCol = mat.diffuseCol.maxColorComponent < mat.emissionCol.maxColorComponent * mat.emissionStrength;
		Color displayCol = displayEmissiveCol ? mat.emissionCol * mat.emissionStrength : mat.diffuseCol;
		displayCol.a = 0.5f;
		meshRenderer.sharedMaterial.color = displayCol;
	}
}