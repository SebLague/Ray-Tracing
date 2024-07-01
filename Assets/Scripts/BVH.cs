using UnityEngine;
using System;
using Unity.Mathematics;
using System.Text;

public class BVH
{
    public readonly NodeList allNodes;
    public readonly Triangle[] allTris;
    public BuildStats stats;

    public Triangle[] GetTriangles() => allTris;

    readonly BVHTriangle[] AllTriangles;

    public BVH(Vector3[] verts, int[] indices, Vector3[] normals)
    {
        // Start recording stats
        var sw = System.Diagnostics.Stopwatch.StartNew();
        stats = new();

        // Construct BVH
        allNodes = new();
        AllTriangles = new BVHTriangle[indices.Length / 3];
        BoundingBox bounds = new BoundingBox();

        for (int i = 0; i < indices.Length; i += 3)
        {
            float3 a = verts[indices[i + 0]];
            float3 b = verts[indices[i + 1]];
            float3 c = verts[indices[i + 2]];
            float3 centre = (a + b + c) / 3;
            float3 max = math.max(math.max(a, b), c);
            float3 min = math.min(math.min(a, b), c);
            AllTriangles[i / 3] = new BVHTriangle(centre, min, max, i);
            bounds.GrowToInclude(min, max);
        }

        allNodes.Add(new Node(bounds));
        Split(0, verts, 0, AllTriangles.Length);

        allTris = new Triangle[AllTriangles.Length];
        for (int i = 0; i < AllTriangles.Length; i++)
        {
            BVHTriangle buildTri = AllTriangles[i];
            Vector3 a = verts[indices[buildTri.Index + 0]];
            Vector3 b = verts[indices[buildTri.Index + 1]];
            Vector3 c = verts[indices[buildTri.Index + 2]];
            Vector3 norm_a = normals[indices[buildTri.Index + 0]];
            Vector3 norm_b = normals[indices[buildTri.Index + 1]];
            Vector3 norm_c = normals[indices[buildTri.Index + 2]];
            allTris[i] = new Triangle(a, b, c, norm_a, norm_b, norm_c);
        }

        // Finish recording stats
        sw.Stop();
        stats.TimeMs = (int)sw.ElapsedMilliseconds;
    }

    void Split(int parentIndex, Vector3[] verts, int triGlobalStart, int triNum, int depth = 0)
    {
        const int MaxDepth = 32;
        Node parent = allNodes.Nodes[parentIndex];
        Vector3 size = parent.CalculateBoundsSize();
        float parentCost = NodeCost(size, triNum);

        (int splitAxis, float splitPos, float cost) = ChooseSplit(parent, triGlobalStart, triNum);

        if (cost < parentCost && depth < MaxDepth)
        {
            BoundingBox boundsLeft = new();
            BoundingBox boundsRight = new();
            int numOnLeft = 0;

            for (int i = triGlobalStart; i < triGlobalStart + triNum; i++)
            {
                BVHTriangle tri = AllTriangles[i];
                if (tri.Centre[splitAxis] < splitPos)
                {
                    boundsLeft.GrowToInclude(tri.Min, tri.Max);

                    BVHTriangle swap = AllTriangles[triGlobalStart + numOnLeft];
                    AllTriangles[triGlobalStart + numOnLeft] = tri;
                    AllTriangles[i] = swap;
                    numOnLeft++;
                }
                else
                {
                    boundsRight.GrowToInclude(tri.Min, tri.Max);
                }
            }

            int numOnRight = triNum - numOnLeft;
            int triStartLeft = triGlobalStart + 0;
            int triStartRight = triGlobalStart + numOnLeft;

            // Split parent into two children
            int childIndexLeft = allNodes.Add(new(boundsLeft, triStartLeft, 0));
            int childIndexRight = allNodes.Add(new(boundsRight, triStartRight, 0));

            // Update parent
            parent.StartIndex = childIndexLeft;
            allNodes.Nodes[parentIndex] = parent;
            stats.RecordNode(depth, false);

            // Recursively split children
            Split(childIndexLeft, verts, triGlobalStart, numOnLeft, depth + 1);
            Split(childIndexRight, verts, triGlobalStart + numOnLeft, numOnRight, depth + 1);
        }
        else
        {
            // Parent is actually leaf, assign all triangles to it
            parent.StartIndex = triGlobalStart;
            parent.TriangleCount = triNum;
            allNodes.Nodes[parentIndex] = parent;
            stats.RecordNode(depth, true, triNum);
        }
    }

    (int axis, float pos, float cost) ChooseSplit(Node node, int start, int count)
    {
        if (count <= 1) return (0, 0, float.PositiveInfinity);

        float bestSplitPos = 0;
        int bestSplitAxis = 0;
        const int numSplitTests = 5;

        float bestCost = float.MaxValue;

        // Estimate best split pos
        for (int axis = 0; axis < 3; axis++)
        {
            for (int i = 0; i < numSplitTests; i++)
            {
                float splitT = (i + 1) / (numSplitTests + 1f);
                float splitPos = Mathf.Lerp(node.BoundsMin[axis], node.BoundsMax[axis], splitT);
                float cost = EvaluateSplit(axis, splitPos, start, count);
                if (cost < bestCost)
                {
                    bestCost = cost;
                    bestSplitPos = splitPos;
                    bestSplitAxis = axis;
                }
            }
        }

        return (bestSplitAxis, bestSplitPos, bestCost);
    }

    float EvaluateSplit(int splitAxis, float splitPos, int start, int count)
    {
        BoundingBox boundsLeft = new();
        BoundingBox boundsRight = new();
        int numOnLeft = 0;
        int numOnRight = 0;

        for (int i = start; i < start + count; i++)
        {
            BVHTriangle tri = AllTriangles[i];
            if (tri.Centre[splitAxis] < splitPos)
            {
                boundsLeft.GrowToInclude(tri.Min, tri.Max);
                numOnLeft++;
            }
            else
            {
                boundsRight.GrowToInclude(tri.Min, tri.Max);
                numOnRight++;
            }
        }

        float costA = NodeCost(boundsLeft.Size, numOnLeft);
        float costB = NodeCost(boundsRight.Size, numOnRight);
        return costA + costB;
    }

    static float NodeCost(Vector3 size, int numTriangles)
    {
        float halfArea = size.x * size.y + size.x * size.z + size.y * size.z;
        return halfArea * numTriangles;
    }

    public struct Node
    {
        public float3 BoundsMin;
        public float3 BoundsMax;
        // Index of first child (if triangle count is negative) otherwise index of first triangle
        public int StartIndex;
        public int TriangleCount;

        public Node(BoundingBox bounds) : this()
        {
            BoundsMin = bounds.Min;
            BoundsMax = bounds.Max;
            StartIndex = -1;
            TriangleCount = -1;
        }

        public Node(BoundingBox bounds, int startIndex, int triCount)
        {
            BoundsMin = bounds.Min;
            BoundsMax = bounds.Max;
            StartIndex = startIndex;
            TriangleCount = triCount;
        }

        public Vector3 CalculateBoundsSize() => BoundsMax - BoundsMin;
        public Vector3 CalculateBoundsCentre() => (BoundsMin + BoundsMax) / 2;
    }

    public struct BoundingBox
    {
        public float3 Min;
        public float3 Max;
        public float3 Centre => (Min + Max) / 2;
        public float3 Size => Max - Min;
        bool hasPoint;

        public void GrowToInclude(float3 min, float3 max)
        {
            if (hasPoint)
            {
                Min.x = min.x < Min.x ? min.x : Min.x;
                Min.y = min.y < Min.y ? min.y : Min.y;
                Min.z = min.z < Min.z ? min.z : Min.z;
                Max.x = max.x > Max.x ? max.x : Max.x;
                Max.y = max.y > Max.y ? max.y : Max.y;
                Max.z = max.z > Max.z ? max.z : Max.z;
            }
            else
            {
                hasPoint = true;
                Min = min;
                Max = max;
            }
        }
    }


    public readonly struct BVHTriangle
    {
        public readonly float3 Centre;
        public readonly float3 Min;
        public readonly float3 Max;
        public readonly int Index;

        public BVHTriangle(float3 centre, float3 min, float3 max, int index)
        {
            Centre = centre;
            Min = min;
            Max = max;
            Index = index;
        }
    }

    public Node[] GetNodes() => allNodes.Nodes.AsSpan(0, allNodes.NodeCount).ToArray();


    public class NodeList
    {
        public Node[] Nodes = new Node[256];
        int Index;

        public int Add(Node node)
        {
            if (Index >= Nodes.Length)
            {
                Array.Resize(ref Nodes, Nodes.Length * 2);
            }

            int nodeIndex = Index;
            Nodes[Index++] = node;
            return nodeIndex;
        }

        public int NodeCount => Index;

    }

    [System.Serializable]
    public class BuildStats
    {
        public int TimeMs;
        public int TriangleCount;
        public int TotalNodeCount;
        public int LeafNodeCount;

        public int LeafDepthMax;
        public int LeafDepthMin = int.MaxValue;
        public int LeafDepthSum;

        public int LeafMaxTriCount;
        public int LeafMinTriCount = int.MaxValue;

        public void RecordNode(int depth, bool isLeaf, int triCount = 0)
        {
            TotalNodeCount++;

            if (isLeaf)
            {
                LeafNodeCount++;
                LeafDepthSum += depth;
                LeafDepthMin = Mathf.Min(LeafDepthMin, depth);
                LeafDepthMax = Mathf.Max(LeafDepthMax, depth);
                TriangleCount += triCount;

                LeafMaxTriCount = Mathf.Max(LeafMaxTriCount, triCount);
                LeafMinTriCount = Mathf.Min(LeafMinTriCount, triCount);
            }
        }


        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"Time (BVH): {TimeMs} ms");
            sb.AppendLine($"Triangles: {TriangleCount}");
            sb.AppendLine($"Node Count: {TotalNodeCount}");
            sb.AppendLine($"Leaf Count: {LeafNodeCount}");
            sb.AppendLine($"Leaf Depth:");
            sb.AppendLine($" - Min: {LeafDepthMin}");
            sb.AppendLine($" - Max: {LeafDepthMax}");
            sb.AppendLine($" - Mean: {LeafDepthSum / (float)LeafNodeCount:0.####}");
            sb.AppendLine($"Leaf Tris:");
            sb.AppendLine($" - Min: {LeafMinTriCount}");
            sb.AppendLine($" - Max: {LeafMaxTriCount}");
            sb.AppendLine($" - Mean: {TriangleCount / (float)LeafNodeCount:0.####}");

            return sb.ToString();
        }
    }
}
