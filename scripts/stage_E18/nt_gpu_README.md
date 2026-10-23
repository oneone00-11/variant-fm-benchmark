# NT 在原 GPU 环境里打分（E18.2）

这个包在一台租来的 GPU 上，用 atlas 当初算 `nt` 列的同一套环境，给 39 个锚点和 6,667 个外部变异打 Nucleotide Transformer 分数。锚点与外部变异在同一次运行里打。锚点 39/39 与分析集逐位相同，外部变异的分数才算数；不是 39/39，NT 就不进外部检验，不设容差。

包里没有任何密钥，也不需要。

## 要什么样的机器

- **GPU：NVIDIA RTX A6000**（Ampere 架构，48 GB）。atlas 的运行记录是 RunPod 上的 A6000（2026-08-01）。不是 RTX 6000 Ada，也不是 A100。脚本会核对显卡名称（`NVIDIA RTX A6000`）和计算能力（8.6），不符就拒绝运行。
- 驱动要支持 CUDA 12.8。驱动版本当初没有记录，这次照实记进 `out/environment.json`。
- 能上网：要从 Hugging Face 下载模型（约 2 GB，不需要账号或令牌）。

## 步骤

1. 把整个包传到机器上，进入包的目录，先核对文件：

   ```
   sha256sum -c MANIFEST.sha256
   ```

2. 照 atlas 当初的做法，建一个只装记录里那些包的环境（Python 3.12）：

   ```
   conda create -n nt python=3.12 -y
   conda activate nt
   pip install --extra-index-url https://download.pytorch.org/whl/cu128 -r requirements_pod.txt
   ```

   `requirements_pod.txt` 就是 atlas 记录的 `pipfreeze_runpod.txt`，只去掉了 conda 自带的那一行 `packaging @ file://...`。

3. 先试吞吐（E18.2 规定每个打分器都先跑锚点和 50 个变异）：

   ```
   python run_nt_gpu.py --trial 2>&1 | tee run_trial.log
   ```

   脚本先核对：包里的文件、显卡、Python 3.12、torch 2.11.0+cu128（CUDA 12.8）、transformers 4.46.3，以及 `pip freeze --all` 与记录逐个相符（只许多出 pip、wheel、packaging），还有模型的版本（`06615c1660c892fc199840c18123f8385b3542a8`）。任何一项不符就停，什么也不打。

   都相符就打 39 个锚点和前 50 个外部变异，写进 `out_trial/`。屏幕最后几行会写 `anchors: k/39 identical` 和全量的推算用时。试跑的分数不用，只看吞吐；锚点若不是 39/39，可以就此停下，不必再跑全量。

4. 全量运行：

   ```
   python run_nt_gpu.py 2>&1 | tee run.log
   ```

   同样先做全部核对，然后打分：先 39 个锚点，再 6,667 个外部变异，同一次运行。按 atlas 的记录，A6000 上每秒约 10–25 个，几分钟到十几分钟。

5. 把 `out/`、`out_trial/` 两个目录和 `run.log`、`run_trial.log` 拿回来。`out/` 里是：
   - `results.tsv`：每个变异的分数（Python 的 repr），锚点那几行带"是否与分析集逐位相同"；
   - `environment.json`：nvidia-smi 输出、驱动、CUDA、cuDNN、显卡名称与计算能力、Python、pip freeze、模型版本、起止时间；
   - `MANIFEST.sha256`：上面这些文件的 sha256。

   屏幕最后一行会写 `anchors: k/39 identical`。拿回来以后，本机会再核对一次校验和与锚点。

## 照原运行、没有试出来的设置

原运行有记录的，一律照记录：显卡、torch、transformers、CUDA、全部包的版本。没有记录的，用 atlas `models/nt/score.py` 的默认值：

- float32；
- 一次只送一个变异；
- 6,000 nt 的窗口；
- 模型仓库默认分支的版本。这个分支自 2025-10-22 起就指向 06615c16，早于原运行。

模型要靠仓库自带的代码才能载入。transformers 只有在交互式回答 "y"，或传入 `trust_remote_code=True` 时，才会执行这段代码。原运行是在 tmux 里交互运行的，只能是答了 "y"，所以这里传 `trust_remote_code=True`。分词器照 score.py，不传任何参数。
