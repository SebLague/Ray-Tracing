Shader "Hidden/Accumulate"
{
	Properties
	{
		_MainTex ("Texture", 2D) = "white" {}
	}
	SubShader
	{
		Cull Off ZWrite Off ZTest Always

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag

			#include "UnityCG.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float2 uv : TEXCOORD0;
				float4 vertex : SV_POSITION;
			};

			v2f vert (appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = v.uv;
				return o;
			}

			sampler2D _MainTex;
			sampler2D _PrevFrame;
			int _Frame;

			float4 frag (v2f i) : SV_Target
			{
				float4 col = tex2D(_MainTex, i.uv);
				float4 colPrev = tex2D(_PrevFrame, i.uv);

				float weight = 1.0 / (_Frame + 1);
				// Combine prev frame with current frame. Weight the contributions to result in an average over all frames.
				float4 accumulatedCol = saturate(colPrev * (1 - weight) + col * weight);
				
				return accumulatedCol;

			}
			ENDCG
		}
	}
}
