#include <benchmark/benchmark.h>
#include "openfhe.h"

using namespace lbcrypto;

// static std::unique_ptr<CCParams<CryptoContextBFVRNS>> parameters;
static std::unique_ptr<CryptoContext<DCRTPoly>> globalCryptoContext;
static std::unique_ptr<KeyPair<DCRTPoly>> globalKeyPair;

//global plaintext
static std::unique_ptr<Plaintext> globalPlaintext1;
static std::unique_ptr<Plaintext> globalPlaintext2;
//global ciphertexts
static std::unique_ptr<Ciphertext<DCRTPoly>> globalCiphertext1;
static std::unique_ptr<Ciphertext<DCRTPoly>> globalCiphertext2;

static int global_slot_count = 0;

static void BFV_Setup(benchmark::State& state) {
    for (auto _ : state) {
        if (!globalCryptoContext) {
            // Your existing setup code...
            CCParams<CryptoContextBFVRNS> parameters;
            parameters.SetSecurityLevel(HEStd_128_classic);
            parameters.SetPlaintextModulus(65537);
            parameters.SetMultiplicativeDepth(2);
            parameters.SetMaxRelinSkDeg(2);
            parameters.SetScalingModSize(20);

            // Store the setup result globally
            globalCryptoContext = std::make_unique<CryptoContext<DCRTPoly>>(GenCryptoContext(parameters));
            globalCryptoContext->get()->Enable(PKE);
            globalCryptoContext->get()->Enable(KEYSWITCH);
            globalCryptoContext->get()->Enable(LEVELEDSHE);
            globalCryptoContext->get()->Enable(ADVANCEDSHE);

            // std::cout << "\np = " << globalCryptoContext->get()->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
            // std::cout << "n = " << globalCryptoContext->get()->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
            //         << std::endl;
            // std::cout << "log2 q = "
            //         << log2(globalCryptoContext->get()->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
            //         << std::endl;

            globalKeyPair = std::make_unique<KeyPair<DCRTPoly>>(globalCryptoContext->get()->KeyGen());    
            globalCryptoContext->get()->EvalMultKeysGen(globalKeyPair->secretKey);

            long slot_count = globalCryptoContext->get()->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2;
            global_slot_count = slot_count;
            // long log_poly_modulus_size = (long)log2(slot_count);
            // std::cout << "log_slot_count = " << log_poly_modulus_size << std::endl;
            // vector<int32_t> rotList(log_poly_modulus_size);


            std::vector<int64_t> plain_vec(slot_count, 0);
            std::vector<int64_t> plain_vec2(slot_count, 0);
            for (int i = 0; i < slot_count; i++) {
                plain_vec[i] = rand() % 2;
                plain_vec2[i] = rand() % 2;
            }

            // encode global Plaintext
            globalPlaintext1 = std::make_unique<Plaintext>(globalCryptoContext->get()->MakeCoefPackedPlaintext(plain_vec));
            globalPlaintext2 = std::make_unique<Plaintext>(globalCryptoContext->get()->MakeCoefPackedPlaintext(plain_vec2));

            globalCiphertext1 = std::make_unique<Ciphertext<DCRTPoly>>(globalCryptoContext->get()->Encrypt(globalKeyPair->publicKey, *globalPlaintext1));
            globalCiphertext2 = std::make_unique<Ciphertext<DCRTPoly>>(globalCryptoContext->get()->Encrypt(globalKeyPair->publicKey, *globalPlaintext2));
        }
    }
}

static void BFV_Enc(benchmark::State& state) {
    std::vector<int64_t> plain_vec(global_slot_count, 0);
    for (int i = 0; i < global_slot_count; i++) {
        plain_vec[i] = rand() % 2;
    }
    Plaintext plaintext = globalCryptoContext->get()->MakeCoefPackedPlaintext(plain_vec);
    for (auto _ : state) {
        // Benchmark this operation
        globalCryptoContext->get()->Encrypt(globalKeyPair->publicKey, plaintext);
    }
}

static void BFV_Dec(benchmark::State& state) {
    std::vector<int64_t> plain_vec(global_slot_count, 0);
    for (int i = 0; i < global_slot_count; i++) {
        plain_vec[i] = rand() % 2;
    }
    Plaintext plaintext = globalCryptoContext->get()->MakeCoefPackedPlaintext(plain_vec);
    Ciphertext<DCRTPoly> ctxt;
    ctxt = globalCryptoContext->get()->Encrypt(globalKeyPair->publicKey, plaintext);

    Plaintext plain_out;
    for (auto _ : state) {
        // Benchmark this operation
        globalCryptoContext->get()->Decrypt(globalKeyPair->secretKey, ctxt, &plain_out);
    }
}

static void BFV_add(benchmark::State& state) {
    for (auto _ : state) {
        // Benchmark this operation
        globalCryptoContext->get()->EvalAdd(*globalCiphertext1, *globalCiphertext2);
    }

}

static void BFV_MultPlain(benchmark::State& state) {
    for (auto _ : state) {
        // Benchmark this operation
        globalCryptoContext->get()->EvalMult(*globalCiphertext1, *globalPlaintext1);
    }
}

static void BFV_Mult(benchmark::State& state) {
    for (auto _ : state) {
        // Benchmark this operation
        globalCryptoContext->get()->EvalMult(*globalCiphertext1, *globalCiphertext2);
    }
}




int main(int argc, char** argv) {
  // Run the benchmark
  benchmark::RegisterBenchmark("BFV_Setup", BFV_Setup);
  benchmark::RegisterBenchmark("BFV_Enc", BFV_Enc);
  benchmark::RegisterBenchmark("BFV_Dec", BFV_Dec);
benchmark::RegisterBenchmark("BFV_add", BFV_add);
benchmark::RegisterBenchmark("BFV_MultPlain", BFV_MultPlain);
benchmark::RegisterBenchmark("BFV_Mult", BFV_Mult);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
// // Register the function as a benchmark
// BENCHMARK(BM_SomeFunction);
// // Run the benchmark
// BENCHMARK_MAIN();