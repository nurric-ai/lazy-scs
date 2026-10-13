"""
Is "fast matrix multiplication" worth it on a real GPU?

The omega <= 9/4 result (openai/math #107) is a non-constructive existence proof,
so there is no algorithm to code. This script tests the one fast-matmul family
that *is* implementable -- Strassen (omega = 2.807) -- which is the best case
for the whole idea: same recursion trick, far smaller hidden constants than
any laser-method / 2.37 / 2.25 scheme.

Three tests:
  1. Speed:     Strassen (various recursion depths) vs torch.matmul
  2. Accuracy:  error vs an fp64 reference, in fp32 and fp16
  3. Training:  a small MLP trained with Strassen linear layers vs normal ones

Run:  python fastmm_bench.py
"""
import math
import time

import torch
import torch.nn as nn
import torch.nn.functional as F

dev = "cuda" if torch.cuda.is_available() else "cpu"
torch.backends.cuda.matmul.allow_tf32 = False  # fair fp32 comparison
torch.manual_seed(0)


# --------------------------------------------------------------------------
# Strassen: 7 sub-multiplications instead of 8 per level. Works for any
# (M,K)x(K,N) with even dims at each level; recursion stops at `depth`.
# --------------------------------------------------------------------------
def strassen(A, B, depth):
    M, K = A.shape[-2:]
    N = B.shape[-1]
    if depth == 0 or M % 2 or K % 2 or N % 2:
        return A @ B
    m, k, n = M // 2, K // 2, N // 2
    A11, A12, A21, A22 = A[..., :m, :k], A[..., :m, k:], A[..., m:, :k], A[..., m:, k:]
    B11, B12, B21, B22 = B[..., :k, :n], B[..., :k, n:], B[..., k:, :n], B[..., k:, n:]
    d = depth - 1
    P1 = strassen(A11 + A22, B11 + B22, d)
    P2 = strassen(A21 + A22, B11, d)
    P3 = strassen(A11, B12 - B22, d)
    P4 = strassen(A22, B21 - B11, d)
    P5 = strassen(A11 + A12, B22, d)
    P6 = strassen(A21 - A11, B11 + B12, d)
    P7 = strassen(A12 - A22, B21 + B22, d)
    C11 = P1 + P4 - P5 + P7
    C12 = P3 + P5
    C21 = P2 + P4
    C22 = P1 - P2 + P3 + P6
    top = torch.cat([C11, C12], dim=-1)
    bot = torch.cat([C21, C22], dim=-1)
    return torch.cat([top, bot], dim=-2)


def timeit(fn, reps=5):
    fn()  # warmup
    if dev == "cuda":
        torch.cuda.synchronize()
    t = time.perf_counter()
    for _ in range(reps):
        fn()
    if dev == "cuda":
        torch.cuda.synchronize()
    return (time.perf_counter() - t) / reps


# --------------------------------------------------------------------------
# 0. The theory: why the paper is exciting on paper
# --------------------------------------------------------------------------
def theory_table():
    print("\n=== 0. Multiplication counts if hidden constants were 1 (they are NOT) ===")
    print(f"{'n':>7} {'n^3':>12} {'Strassen 2.807':>16} {'2.371':>12} {'9/4':>12}  speedup n^3/n^2.25")
    for n in [1024, 4096, 16384]:
        print(f"{n:>7} {n**3:>12.2e} {n**2.807:>16.2e} {n**2.371:>12.2e} {n**2.25:>12.2e}  {n**0.75:>8.0f}x")
    print("-> looks like 500x+ for free. The rest of this script measures reality.")


# --------------------------------------------------------------------------
# 1 + 2. Speed and accuracy
# --------------------------------------------------------------------------
def speed_and_accuracy():
    print(f"\n=== 1/2. Speed + accuracy on {torch.cuda.get_device_name() if dev=='cuda' else 'CPU'} ===")
    sizes = [1024, 2048, 4096, 8192]
    for dtype in [torch.float32, torch.float16]:
        print(f"\n--- {str(dtype).replace('torch.', '')} ---")
        print(f"{'n':>6} {'depth':>5} {'time ms':>9} {'vs matmul':>10} {'rel err':>10}")
        for n in sizes:
            A64 = torch.randn(n, n, device=dev, dtype=torch.float64)
            B64 = torch.randn(n, n, device=dev, dtype=torch.float64)
            ref = A64 @ B64
            A, B = A64.to(dtype), B64.to(dtype)
            del A64, B64
            base_t = None
            for depth in [0, 1, 2, 3]:
                if n >> depth < 256:
                    continue
                f = lambda: strassen(A, B, depth)
                t = timeit(f, reps=3 if n >= 8192 else 5)
                err = ((f().double() - ref).norm() / ref.norm()).item()
                if depth == 0:
                    base_t = t
                tag = "matmul" if depth == 0 else f"{base_t / t:.2f}x"
                print(f"{n:>6} {depth:>5} {t*1e3:>9.2f} {tag:>10} {err:>10.2e}")
            del A, B, ref
            if dev == "cuda":
                torch.cuda.empty_cache()


# --------------------------------------------------------------------------
# 3. Does it hurt training?  Small MLP on a synthetic teacher task.
# --------------------------------------------------------------------------
class StrassenLinear(nn.Linear):
    def __init__(self, i, o, depth):
        super().__init__(i, o)
        self.depth = depth

    def forward(self, x):
        return strassen(x, self.weight.t(), self.depth) + self.bias


def make_mlp(width, depth_s):
    L = (lambda i, o: StrassenLinear(i, o, depth_s)) if depth_s else nn.Linear
    return nn.Sequential(L(512, width), nn.GELU(), L(width, width), nn.GELU(), L(width, 16)).to(dev)


def training_test(steps=300, width=4096, batch=4096):
    print(f"\n=== 3. Training: MLP 512->{width}->{width}->16, batch {batch}, fp16 autocast, {steps} steps ===")
    g = torch.Generator(device=dev).manual_seed(1)
    teacher = torch.randn(512, 16, device=dev, generator=g) / math.sqrt(512)
    for depth_s in [0, 1, 2]:
        torch.manual_seed(0)
        model = make_mlp(width, depth_s)
        opt = torch.optim.AdamW(model.parameters(), lr=1e-3)
        gx = torch.Generator(device=dev).manual_seed(2)

        def step():
            x = torch.randn(batch, 512, device=dev, generator=gx)
            y = torch.sin(x @ teacher)
            with torch.autocast(dev, dtype=torch.float16):
                loss = F.mse_loss(model(x).float(), y)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            return loss

        step()  # warmup
        if dev == "cuda":
            torch.cuda.synchronize()
        t = time.perf_counter()
        for _ in range(steps):
            loss = step()
        if dev == "cuda":
            torch.cuda.synchronize()
        dt = (time.perf_counter() - t) / steps
        name = "nn.Linear" if depth_s == 0 else f"Strassen d={depth_s}"
        print(f"{name:>14}: {dt*1e3:7.2f} ms/step   final loss {loss.item():.5f}")


if __name__ == "__main__":
    theory_table()
    speed_and_accuracy()
    training_test()
