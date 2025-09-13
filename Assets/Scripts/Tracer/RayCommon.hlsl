// --- Settings and constants ---
static const float PI = 3.1415;

// Raytracing Settings
int MaxBounceCount;
int NumRaysPerPixel;
int Frame;
int renderSeed;

// Camera settings
float DefocusStrength;
float DivergeStrength;
float3 ViewParams;
float4x4 CamLocalToWorldMatrix;

// Sky settings
int UseSky;
float3 SunColour;
float SunFocus = 500;
float SunIntensity = 10;
float3 dirToSun;

// Debug settings
int visMode;
float debugVisScale;
float4 debugParams;

// Constants
static const int MATERIAL_CHECKERED = 1;
static const int MATERIAL_GLASS = 2;
static const float EPSILON_GLASS = 0.01;

// ---- Structures ----

struct Ray
{
    // Current position and direction of the ray
    // (inverse dir cached to speed up collision tests)
    float3 pos;
    float3 dir;
    float3 invDir;
    // Proportion of rgb light that can be transmitted along
    // this ray's path (if it reaches a light source)
    float3 transmittance;
    // Number of times the path has collided with an object
    int bounceCount;
};

struct Triangle
{
    float3 posA, posB, posC;
    float3 normA, normB, normC;
};

struct TriangleHitInfo
{
    bool didHit;
    bool isBackface;
    float dst;
    float3 hitPoint;
    float3 normal;
};

struct RayTracingMaterial
{
    float4 diffuseCol;
    float4 emissionCol;
    float4 specularCol;
    float4 absorption;
    float absorptionStrength;
    float emissionStrength;
    float smoothness;
    float specularProbability;
    float ior;
    int flag;
};

struct Model
{
    int nodeOffset;
    int triOffset;
    float4x4 worldToLocalMatrix;
    float4x4 localToWorldMatrix;
    RayTracingMaterial material;
};

struct BVHNode
{
    float3 boundsMin;
    float3 boundsMax;
    // index refers to triangles if is leaf node (triangleCount > 0)
    // otherwise it is the index of the first child node
    int startIndex;
    int triangleCount;
};

struct ModelHitInfo
{
    bool didHit;
    bool isBackface;
    float3 normal;
    float3 pos;
    float dst;
    RayTracingMaterial material;
};

struct LightResponse
{
    float3 reflectDir;
    float3 refractDir;
    float reflectWeight;
    float refractWeight;
};

// --- Buffers (and their sizes) ---	
StructuredBuffer<Model> ModelInfo;
StructuredBuffer<Triangle> Triangles;
StructuredBuffer<BVHNode> Nodes;
AppendStructuredBuffer<float3> DebugData; //
int triangleCount;
int modelCount;

// ---- RNG Functions ----

// PCG (permuted congruential generator). Thanks to:
// www.pcg-random.org and www.shadertoy.com/view/XlGcRh
uint NextRandom(inout uint state)
{
    state = state * 747796405 + 2891336453;
    uint result = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
    result = (result >> 22) ^ result;
    return result;
}

float RandomValue(inout uint state)
{
    return NextRandom(state) / 4294967295.0; // 2^32 - 1
}

// Random value in normal distribution (with mean=0 and sd=1)
float RandomValueNormalDistribution(inout uint state)
{
    // Thanks to https://stackoverflow.com/a/6178290
    float theta = 2 * 3.1415926 * RandomValue(state);
    float rho = sqrt(-2 * log(RandomValue(state)));
    return rho * cos(theta);
}

// Calculate a random direction
float3 RandomDirection(inout uint state)
{
    // Thanks to https://math.stackexchange.com/a/1585996
    float x = RandomValueNormalDistribution(state);
    float y = RandomValueNormalDistribution(state);
    float z = RandomValueNormalDistribution(state);
    return normalize(float3(x, y, z));
}

float2 RandomPointInCircle(inout uint rngState)
{
    float angle = RandomValue(rngState) * 2 * PI;
    float2 pointOnCircle = float2(cos(angle), sin(angle));
    return pointOnCircle * sqrt(RandomValue(rngState));
}

// Crude sky colour function for background light
float3 GetEnvironmentLight(float3 dir)
{
    if (UseSky == 0) return 0;
    const float3 GroundColour = float3(0.35, 0.3, 0.35);
    const float3 SkyColourHorizon = float3(1, 1, 1);
    const float3 SkyColourZenith = float3(0.08, 0.37, 0.73);


    float skyGradientT = pow(smoothstep(0, 0.4, dir.y), 0.35);
    float groundToSkyT = smoothstep(-0.01, 0, dir.y);
    float3 skyGradient = lerp(SkyColourHorizon, SkyColourZenith, skyGradientT);
    float s = 1000 * 1 / SunFocus;
    float sun = pow(max(0, dot(dir, dirToSun)), s) * SunIntensity;
    // Combine ground, sky, and sun
    float3 composite = lerp(GroundColour, skyGradient, groundToSkyT) + sun * SunColour * (groundToSkyT >= 1);
    return composite;
}

// --- Ray Intersection Functions ---

// Ray-triangle intersection (thanks to stackoverflow.com/a/42752998)
TriangleHitInfo RayTriangle(Ray ray, Triangle tri, bool cullBackface)
{
    float3 edgeAB = tri.posB - tri.posA;
    float3 edgeAC = tri.posC - tri.posA;
    float3 triFaceVector = cross(edgeAB, edgeAC);
    float3 vertRayOffset = ray.pos - tri.posA;
    float3 rayOffsetPerp = cross(vertRayOffset, ray.dir);
    float determinant = -dot(ray.dir, triFaceVector);
    float invDet = 1 / determinant;

    // Calculate hit-dst and barycentric coordinates
    float dst = dot(vertRayOffset, triFaceVector) * invDet;
    float u = dot(edgeAC, rayOffsetPerp) * invDet;
    float v = -dot(edgeAB, rayOffsetPerp) * invDet;
    float w = 1 - u - v;

    // Initialize hit info
    TriangleHitInfo hitInfo;
    bool keep = cullBackface ? determinant >= 1E-8 : abs(determinant) >= 1E-8;
    hitInfo.didHit = keep && dst > 0 && u >= 0 && v >= 0 && w >= 0;
    float3 smoothNormal = normalize(tri.normA * w + tri.normB * u + tri.normC * v);
    hitInfo.normal = smoothNormal * sign(determinant);
    // hitInfo.normal = normalize(triFaceVector);
    hitInfo.isBackface = determinant < 0;
    hitInfo.hitPoint = ray.pos + ray.dir * dst;
    hitInfo.dst = dst;
    return hitInfo;
}


// Thanks to https://tavianator.com/2011/ray_box.html
float RayBoundingBoxDst(Ray ray, float3 boxMin, float3 boxMax)
{
    float3 tMin = (boxMin - ray.pos) * ray.invDir;
    float3 tMax = (boxMax - ray.pos) * ray.invDir;
    float3 t1 = min(tMin, tMax);
    float3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);

    bool hit = tFar >= tNear && tFar > 0;
    float dst = hit ? tNear > 0 ? tNear : 0 : 1.#INF;
    return dst;
};


TriangleHitInfo RayTriangleBVH(inout Ray ray, float rayLength, int nodeOffset, int triOffset, inout int2 stats, bool cullBackface)
{
    TriangleHitInfo result;
    result.dst = rayLength;

    int stack[32];
    int stackCount = 0;
    stack[stackCount++] = nodeOffset + 0;

    while (stackCount > 0)
    {
        BVHNode node = Nodes[stack[--stackCount]];
        bool isLeaf = node.triangleCount > 0;

        if (isLeaf)
        {
            for (int i = 0; i < node.triangleCount; i++)
            {
                Triangle tri = Triangles[triOffset + node.startIndex + i];
                TriangleHitInfo triHitInfo = RayTriangle(ray, tri, cullBackface);
                stats[0]++; // count triangle intersection tests

                if (triHitInfo.didHit && triHitInfo.dst < result.dst)
                {
                    result = triHitInfo;
                }
            }
        }
        else
        {
            int childIndexA = nodeOffset + node.startIndex + 0;
            int childIndexB = nodeOffset + node.startIndex + 1;
            BVHNode childA = Nodes[childIndexA];
            BVHNode childB = Nodes[childIndexB];

            float dstA = RayBoundingBoxDst(ray, childA.boundsMin, childA.boundsMax);
            float dstB = RayBoundingBoxDst(ray, childB.boundsMin, childB.boundsMax);
            stats[1] += 2; // count bounding box intersection tests

            // We want to look at closest child node first, so push it last
            bool isNearestA = dstA <= dstB;
            float dstNear = isNearestA ? dstA : dstB;
            float dstFar = isNearestA ? dstB : dstA;
            int childIndexNear = isNearestA ? childIndexA : childIndexB;
            int childIndexFar = isNearestA ? childIndexB : childIndexA;

            if (dstFar < result.dst) stack[stackCount++] = childIndexFar;
            if (dstNear < result.dst) stack[stackCount++] = childIndexNear;
        }
    }


    return result;
}

ModelHitInfo RaySphere(float3 rayPos, float3 rayDir, float3 sphereCentre, float sphereRadius)
{
    ModelHitInfo hitInfo = (ModelHitInfo)0;
    hitInfo.dst = 1.#INF;

    float3 offsetRayOrigin = rayPos - sphereCentre;
    // From the equation: sqrLength(rayOrigin + rayDir * dst) = radius^2
    // Solving for dst results in a quadratic equation with coefficients:
    float a = dot(rayDir, rayDir); // a = 1 (assuming unit vector)
    float b = 2 * dot(offsetRayOrigin, rayDir);
    float c = dot(offsetRayOrigin, offsetRayOrigin) - sphereRadius * sphereRadius;
    // Quadratic discriminant
    float discriminant = b * b - 4 * a * c;

    // No solution when d < 0 (ray misses sphere)
    if (discriminant >= 0)
    {
        float s = sqrt(discriminant);
        // Distance to nearest intersection point (from quadratic formula)
        float dstNear = max(0, (-b - s) / (2 * a));
        float dstFar = (-b + s) / (2 * a);

        // Ignore intersections that occur behind the ray
        if (dstFar >= 0)
        {
            hitInfo.didHit = true;
            bool isInside = dstNear == 0;
            hitInfo.isBackface = isInside;
            hitInfo.dst = isInside ? dstFar : dstNear;

            hitInfo.pos = rayPos + rayDir * hitInfo.dst;
            hitInfo.normal = normalize(hitInfo.pos - sphereCentre) * (isInside ? -1 : 1);
            RayTracingMaterial mat = (RayTracingMaterial)0;
            mat.flag = 3;
            mat.diffuseCol = 1;
            mat.smoothness = 0;
            mat.ior = 1.6;

            hitInfo.material = mat;
        }
    }

    return hitInfo;
}


ModelHitInfo CalculateRayCollision(Ray worldRay, bool forceDontCullBack)
{
    ModelHitInfo result;
    result.dst = 1.#INF;
    int2 stats;

    //result = RaySphere(worldRay.origin, worldRay.dir, float3(0, 1.8, 0), 1);

    Ray localRay;
    localRay.transmittance = 0;
    localRay.bounceCount = 0;

    for (int i = 0; i < modelCount; i++)
    {
        Model model = ModelInfo[i];
        // Transform ray into model's local coordinate space
        localRay.pos = mul(model.worldToLocalMatrix, float4(worldRay.pos, 1)).xyz;
        localRay.dir = mul(model.worldToLocalMatrix, float4(worldRay.dir, 0)).xyz;
        localRay.invDir = 1 / localRay.dir; //

        bool cullBackface = model.material.flag != MATERIAL_GLASS;
        if (forceDontCullBack) cullBackface = false;
        // cullBackface = true;
        // Traverse bvh to find closest triangle intersection with current model
        TriangleHitInfo hit = RayTriangleBVH(localRay, result.dst, model.nodeOffset, model.triOffset, stats, cullBackface);

        // Record closest hit
        if (hit.dst < result.dst)
        {
            result.didHit = true;
            result.isBackface = hit.isBackface;
            result.dst = hit.dst;
            result.normal = normalize(mul(model.localToWorldMatrix, float4(hit.normal, 0))).xyz;
            result.pos = worldRay.pos + worldRay.dir * hit.dst;
            result.material = model.material;
        }
    }

    return result;
}

float2 mod2(float2 x, float2 y)
{
    return x - y * floor(x / y);
}

// Calculate the proportion of light that is reflected at the boundary between two media (via the fresnel equations)
// Note: the amount of light refracted can be calculated as 1 minus this value
float CalculateReflectance(float3 inDir, float3 normal, float iorA, float iorB)
{
    float refractRatio = iorA / iorB;
    float cosAngleIn = -dot(inDir, normal);
    float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
    if (sinSqrAngleOfRefraction >= 1) return 1; // Ray is fully reflected, no refraction occurs

    float cosAngleOfRefraction = sqrt(1 - sinSqrAngleOfRefraction);
    float denominatorPerpendicular = iorA * cosAngleIn + iorB * cosAngleOfRefraction;
    float denominatorParallel = iorA * cosAngleIn + iorB * cosAngleOfRefraction;

    if (min(denominatorPerpendicular, denominatorParallel) < 1E-8) return 1;

    // Perpendicular polarization
    float rPerpendicular = (iorA * cosAngleIn - iorB * cosAngleOfRefraction) / denominatorPerpendicular;
    rPerpendicular *= rPerpendicular;
    // Parallel polarization
    float rParallel = (iorB * cosAngleIn - iorA * cosAngleOfRefraction) / denominatorParallel;
    rParallel *= rParallel;

    // Return the average of the perpendicular and parallel polarizations
    return (rPerpendicular + rParallel) / 2;
}


float3 Refract(float3 inDir, float3 normal, float iorA, float iorB)
{
    float refractRatio = iorA / iorB;
    float cosAngleIn = -dot(inDir, normal);
    float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
    if (sinSqrAngleOfRefraction > 1) return 0; // Ray is fully reflected, no refraction occurs

    float3 refractDir = refractRatio * inDir + (refractRatio * cosAngleIn - sqrt(1 - sinSqrAngleOfRefraction)) * normal;
    return refractDir;
}

float3 Reflect(float3 inDir, float3 normal)
{
    return inDir - 2 * dot(inDir, normal) * normal;
}

LightResponse CalculateReflectionAndRefraction(float3 inDir, float3 normal, float iorA, float iorB)
{
    LightResponse result;

    // Calculate the two directions that light can take
    result.reflectDir = Reflect(inDir, normal);
    result.refractDir = Refract(inDir, normal, iorA, iorB);

    // Calculate the proportion of light [0, 1] that takes each path
    result.reflectWeight = CalculateReflectance(inDir, normal, iorA, iorB);
    result.refractWeight = 1 - result.reflectWeight;

    return result;
}

Ray CreateRay(float3 origin, float3 dir, float3 transmittance, int bounceIndex)
{
    Ray ray;
    ray.pos = origin;
    ray.dir = dir;
    ray.invDir = 1 / dir;
    ray.transmittance = transmittance;
    ray.bounceCount = bounceIndex;
    return ray;
}

float3 GetMaterialColour(RayTracingMaterial mat, float3 pos, float3 normal, bool isSpecularBounce)
{
    float3 col = mat.diffuseCol.rgb;

    if (mat.flag == MATERIAL_CHECKERED)
    {
        float2 checkerPoint = pos.xz;
        if (abs(normal.x) > abs(normal.y)) checkerPoint = pos.zy;
        if (abs(normal.z) > max(abs(normal.x), abs(normal.y))) checkerPoint = pos.xy;
        
        checkerPoint *= 1.5;
        float2 c = mod2(floor(checkerPoint), 2.0);
        col = c.x == c.y ? col : mat.emissionCol.rgb;
    }

    return lerp(col, mat.specularCol.rgb, isSpecularBounce);
}

static const float epsilon = 0.001;


// Calculate a random direction
float3 RandomHemisphereDirection(float3 normal, inout uint state)
{
    float3 dir = RandomDirection(state);
    return dir * sign(dot(normal, dir));
}


float3 Trace(Ray initialRay, inout uint rngState)
{
    float3 totalLight = 0;
    Ray ray = initialRay;

    // Bounce the ray around the world to gather light
    for (int i = ray.bounceCount; i <= MaxBounceCount; i++)
    {
        ModelHitInfo hit = CalculateRayCollision(ray, false);
        if (!hit.didHit)
        {
            if (UseSky)
            {
                totalLight += ray.transmittance * GetEnvironmentLight(ray.dir);
            }
            break;
        }

        RayTracingMaterial material = hit.material;

        if (material.flag == MATERIAL_GLASS) // Glass-like material
        {
            // Absorb some amount of light as it travels through the object
            if (hit.isBackface) ray.transmittance *= exp(-hit.dst * material.absorption.rgb * material.absorptionStrength);

            float iorCurrent = hit.isBackface ? material.ior : 1;
            float iorNext = hit.isBackface ? 1 : material.ior;
            LightResponse lr = CalculateReflectionAndRefraction(ray.dir, hit.normal, iorCurrent, iorNext);

            // Calculate random direction in hemisphere around surface normal (cosine-weighted)
            float3 diffuseDir = normalize(hit.normal + RandomDirection(rngState));
            // Randomize the reflect/refract directions based on smoothness for a frosted effect
            lr.reflectDir = normalize(lerp(diffuseDir, lr.reflectDir, material.specularProbability));
            lr.refractDir = normalize(lerp(-diffuseDir, lr.refractDir, material.smoothness));

            // Choose between reflection and refraction probabilistically based on proportion of light going each way
            bool followReflection = RandomValue(rngState) <= lr.reflectWeight;
            ray.dir = followReflection ? lr.reflectDir : lr.refractDir;
            ray.pos = hit.pos + epsilon * hit.normal * sign(dot(hit.normal, ray.dir));
        }
        else
        {
            bool isSpecularBounce = material.specularProbability >= RandomValue(rngState);

            // Redirect ray based on collision info
            ray.pos = hit.pos + (hit.normal * epsilon);
            float3 diffuseDir = normalize(hit.normal + RandomDirection(rngState));
            float3 specularDir = reflect(ray.dir, hit.normal);
            ray.dir = normalize(lerp(diffuseDir, specularDir, material.smoothness * isSpecularBounce));

            // Update light info
            float3 emittedLight = material.emissionCol.rgb * material.emissionStrength;
            totalLight += emittedLight * ray.transmittance;
            ray.transmittance *= GetMaterialColour(material, hit.pos, hit.normal, isSpecularBounce);
        }

        // Randomly early-exit paths, with probability based on how little light can be transmitted along it
        float p = max(ray.transmittance.r, max(ray.transmittance.g, ray.transmittance.b));
        if (RandomValue(rngState) >= p) break;
        ray.transmittance *= 1 / p; // scale by inverse probability so result averages out over many iterations
    }

    return totalLight;
}


float3 RayTrace(float2 uv, uint2 numPixels)
{
    float3 camOrigin = mul(CamLocalToWorldMatrix, float4(0, 0, 0, 1)).xyz;

    // Create seed for random number generator
    uint2 pixelCoord = uv * numPixels;
    uint pixelIndex = pixelCoord.y * numPixels.x + pixelCoord.x;
    uint rngState = pixelIndex + Frame * 719393 + renderSeed; //

    // Calculate focus point
    float3 focusPointLocal = float3(uv - 0.5, 1) * ViewParams;
    float3 focusPoint = mul(CamLocalToWorldMatrix, float4(focusPointLocal, 1)).xyz;
    float3 camRight = CamLocalToWorldMatrix._m00_m10_m20;
    float3 camUp = CamLocalToWorldMatrix._m01_m11_m21;

    // Trace multiple rays and average together
    float3 totalIncomingLight = 0;

    for (int rayIndex = 0; rayIndex < NumRaysPerPixel; rayIndex++)
    {
        // -- Calculate ray origin and direction --
        // Jitter the starting point of the ray. This allows for a depth of field effect.
        float2 defocusJitter = RandomPointInCircle(rngState) * DefocusStrength / numPixels.x;
        float3 rayOrigin = camOrigin + camRight * defocusJitter.x + camUp * defocusJitter.y;

        // Jitter the focus point when calculating the ray direction to allow for blurring the image
        // (at low strengths, this can be used for anti-aliasing)
        float2 jitter = RandomPointInCircle(rngState) * DivergeStrength / numPixels.x;
        float3 jitteredFocusPoint = focusPoint + camRight * jitter.x + camUp * jitter.y;
        float3 rayDir = normalize(jitteredFocusPoint - rayOrigin);

        Ray ray = CreateRay(rayOrigin, rayDir, 1, 0);

        totalIncomingLight += Trace(ray, rngState);
    }

    return totalIncomingLight / NumRaysPerPixel;
}
