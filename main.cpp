#include <Novice.h>
#include "MatrixFunction.h"
#include <chrono>

const char kWindowTitle[] = "GC2B_14_ハラ_アイ";

static const int kRowHeight = 20;
static const int kColumnWidth = 60;

const int kWindowWidth = 1280;
const int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	Vector3 v1{ 1.2f, -3.9f, 2.5f };
	Vector3 v2{ 2.8f, 0.4f, -1.3f };
	Vector3 cameraPosition{ 0.0f,0.0f,-5.0f };
	Vector3 rotate{ 0.0f,0.0f,0.0f };
	Vector3 translate{};

	Vector3 cross = Cross(v1, v2);

	Vector3 vertices[3] =
	{
		{0.0f, 0.5f, 0.0f},
		{-0.5f, -0.5f, 0.0f},
		{0.5f, -0.5f, 0.0f}
	};


	auto start = std::chrono::high_resolution_clock::now();

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		//移動の更新処理
		if (keys[DIK_W])
		{
			translate.z += 0.01f;
		}

		if (keys[DIK_S])
		{
			translate.z -= 0.01f;
		}

		if (keys[DIK_A])
		{
			translate.x -= 0.01f;
		}


		if (keys[DIK_D])
		{
			translate.x += 0.01f;
		}

		// Y軸回転
		auto now = std::chrono::high_resolution_clock::now();
		float deltaTime = std::chrono::duration<float>(now - start).count();
		rotate.y = deltaTime;

		// 各種行列の計算
		Matrix4x4 worldMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix({ 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 0.0f }, cameraPosition);
		Matrix4x4 viewMatrix = Inverse(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 worldViewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);
		Vector3 screenVertices[3] = {};

		//for (uint32_t i = 0; i < 3; ++i)
		//{
		//	Vector3 ndcVertex = Transform(vertices[i], worldViewProjectionMatrix);
		//	screenVertices[i] = Transform(ndcVertex, viewportMatrix);
		//}

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///


		Novice::DrawTriangle
		(
			int(screenVertices[0].x), int(screenVertices[0].y),
			int(screenVertices[1].x), int(screenVertices[1].y),
			int(screenVertices[2].x), int(screenVertices[2].y),
			RED, kFillModeSolid
		);


		VectorScreenPrintf(0, 0, cross, "Cross");

		MatrixScreenPrintf(0, 30, worldViewProjectionMatrix);
		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}

