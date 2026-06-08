S
n リポジトリ `energy-loss` を作成したい。

目的は、荷電粒子が物質中を通過するときのエネルギー損失、stopping power、range-energy relation を計算する小さなライブラリを作ること。

主な用途は2つある。

1. スペクトロスコピー実験におけるエネルギー損失補正

   * 標的、検出器ガス、窓材、シンチレータ、空気、He bag などを通過するときの平均エネルギー損失を計算する。
   * 粒子、初期運動エネルギー、物質、厚さを与えて、dE/dx、ΔE、出口エネルギーを計算する。
   * 複数 layer からなる material stack も扱えるようにする。

2. 原子核乾板における range → energy 変換

   * 乾板中の短い飛跡、たとえば 10 um 程度の range から粒子の運動エネルギーを求める。
   * この用途では単純な Bethe formula だけに依存せず、range-energy table を読み込んで補間できる設計にする。

重要な設計方針:

* リポジトリ名は `energy-loss`
* Python package 名は `energy_loss`
* Docker は使わない。
* まずは小さく、読みやすく、テストしやすい実装にする。
* 最初のモデルは heavy charged particles 用の basic Bethe formula とする。
* ただし、将来的に density-effect correction, shell correction, Barkas correction, Bloch correction, effective charge, nuclear stopping, NIST/SRIM/ATIMA table comparison を追加できる構造にする。
* SRIM や ATIMA のソースコードは含めない。外部参照テーブルとして比較できる設計にする。
* 単位は明示的に扱う。内部単位を決め、docstring に必ず書く。
* 物理定数や粒子質量は一箇所にまとめる。
* まずは陽子、π±、K±、μ±、電子は粒子定義だけ用意する。ただし v0.1 の Bethe 計算対象は heavy charged particles とし、電子の stopping power は未実装でよい。
* material は最初に H2, LH2, He, air, P10, plastic scintillator, kapton, mylar, aluminum, carbon, nuclear_emulsion を入れる。
* nuclear_emulsion は暫定組成・密度として実装し、README で「実験・乾板ごとの calibration table を使うべき」と明記する。

実装してほしい v0.1:

1. `pyproject.toml`

   * Python 3.11 以上
   * 依存はなるべく少なくする
   * runtime dependency は numpy 程度
   * dev dependency として pytest, ruff を入れる

2. ディレクトリ構造
   energy_loss/
   **init**.py
   constants.py
   particles.py
   materials.py
   units.py
   stopping/
   **init**.py
   bethe.py
   models.py
   transport/
   **init**.py
   layer.py
   propagate.py
   range/
   **init**.py
   csda.py
   table.py
   emulsion/
   **init**.py
   range_energy.py
   tests/
   test_bethe.py
   test_transport.py
   test_range_table.py
   examples/
   spectroscopy_loss.py
   emulsion_range_to_energy.py
   README.md
   AGENTS.md

3. API のイメージ

```python
from energy_loss import Particle, Material, Layer
from energy_loss.transport import propagate

layers = [
    Layer(material="kapton", thickness=50.0, unit="um"),
    Layer(material="P10", thickness=30.0, unit="cm"),
]

result = propagate(
    particle="proton",
    kinetic_energy=200.0,
    layers=layers,
    energy_unit="MeV",
    model="bethe",
)

print(result.exit_energy)
print(result.total_energy_loss)
```

range table の API:

```python
from energy_loss.range import RangeEnergyTable

table = RangeEnergyTable.from_csv(
    "examples/data/emulsion_alpha.csv",
    range_unit="um",
    energy_unit="MeV",
)

energy = table.energy_from_range(10.0)
```

4. Bethe formula

* heavy charged particle 用の基本形を実装する。
* 入力は kinetic energy, particle mass, charge, material Z/A, mean excitation energy I, density。
* 出力は mass stopping power [MeV cm2/g] と linear stopping power [MeV/cm]。
* density effect correction などは v0.1 では optional placeholder でよい。
* 低エネルギーで Bethe formula が破綻しうることを docstring と README に書く。
* 物理的に危険な領域では warning を出す設計にする。たとえば beta が小さすぎる場合など。

5. transport

* finite thickness の layer を通した energy loss を数値的に計算する。
* 厚さを小さい step に分けて、各 step で stopping power を更新する簡単な方法でよい。
* 粒子が途中で止まった場合を扱う。
* result object には initial_energy, exit_energy, total_energy_loss, stopped, layers を含める。

6. range

* CSDA range を numerical integration で計算する関数を作る。
* range → energy は table interpolation を第一にする。
* Bethe formula による inverse range は v0.1 では簡易実装でよい。
* RangeEnergyTable は CSV を読み込み、range_from_energy と energy_from_range を線形補間で返す。

7. tests

* 単位変換のテスト
* Bethe stopping power が正の値を返すテスト
* 厚さゼロなら energy loss がゼロになるテスト
* layer を通すと energy が減るテスト
* RangeEnergyTable の補間テスト
* 粒子が止まるケースのテスト

8. README
   README には以下を書く。

* プロジェクトの目的
* 2つの用途:

  1. spectroscopy experiment の detector/target energy loss correction
  2. nuclear emulsion の range-to-energy conversion
* v0.1 の制限
* Bethe formula は heavy charged particle の比較的高エネルギー領域向けであり、低エネルギー・短飛程では table-based range-energy relation が望ましいこと
* SRIM, ATIMA, NIST STAR などは外部比較対象であり、このリポジトリには含めないこと
* 簡単な使用例

9. AGENTS.md
   AGENTS.md には、このリポジトリで自動 coding agent が守るべき方針を書く。

* 物理式を変更するときは README または docstring に根拠を書く
* 単位を曖昧にしない
* テストなしで物理ロジックを変更しない
* SRIM/ATIMA のコードをコピーしない
* Docker を追加しない
* 大きな依存関係を追加しない
* public API を壊す場合は README を更新する

作業手順:

1. まず現在のディレクトリを確認する。
2. まだファイルがない前提で、上記のリポジトリを初期化する。
3. 実装後、`pytest` と `ruff check` を実行する。
4. 失敗したら修正する。
5. 最後に、作成・変更したファイル一覧、実装した機能、未実装の制限、次にやるべきことを簡潔に報告する。

