# 简介
Mikeno 是 OpenFOAM 的魔改版（分支），面向化工应用。

[简体中文] [English](./README-EN.md)

## 这个分支的目标
1. 更好的支持化工CFD模拟
2. 持续合并上游（基金会版本）的更新
3. 构建系统换血：用cmake全面替代wmake（已完成）
4. 修复编译器优化后意外触发SIGFPE的bug，主要针对单精度。
   - 编译器会用向量化指令加速浮点运算，但有时会产生不应该出现的NAN。这些NAN值不会被使用，却触发SIGFPE。

## 修改内容

### 新工具：
1. `autoCompress`: 自动压缩较大的流场文件（非均匀的场，包括网格）。不论是ascii还是binary都压缩为gzip格式（其他格式不能被OpenFOAM和paraview直接读取）

### 真正的CMake构建系统：
1. 将所有子项目的wmake构建系统替换为Modern cmake。
   - 与现代IDE更加兼容，clangd提供良好的语法提示。
   - 更方便引入外部依赖
   - 更高效的并行编译
   - `find_package(Mikeno CONFIG)`即可导入Mikeno，方便二次开发。
2. 性能优化
   - `Opt`和`Prof`模式编译时开启`-march=native`，允许编译器充分利用SIMD。
3. 支持 `AOCC`
4. 删除一些求解器模块代码中的冗余引用。
   1. 目前清理了 `isothermalFluidSolver`、`fluidSolver`、`multicomponentFluidSolver`、`XiFluidSolver`

### 支持任意高的压力
1. 向 `physicalProperties` 字典中添加 `pOffset` 字段（可选），允许在压速耦合计算时采用表压，更新物性时用绝对压力。**这样较小的压降不会再被浮点计算误差掩盖（特别是单精度模式）**

### 多孔介质传热
1. `porousMediaFluidSolver` 模块允许任意数量的多孔相，可以模拟流-固传热和固-固传热，同时支持热平衡模式和热非平衡模式。

### 更严格的热力学
1. 状态方程：增加 `RedlichKwongGas`，重写 `PengRobinsonGas`
2. 以上的真实气体状态方程都用AspenPlusV14检验过
3. Van der vaals 混合规则
   - `PengRobinsonGas` 支持二元交互系数 `k_ij`
4. 所有热力学模型（比如 `hConst`）具有形如 `ideal_*` 的接口，提供理想气体的物性。
5. 新的混合物模型 `realGasMulticomponentMixture`
   1. 依靠混合气体的混合状态方程计算 `rho`
   2. 从剩余性质计算以下物性的真实性质： `Cp` `Cv` `hs` `ha` `es` `ea`。所有剩余性质都是从混合状态方程计算的。
6. 新的混合物模型 `idealLiquidMulticomponentMixture`
   1. 体积加权平均计算密度（无混合体积）
   2. 质量加权平均计算比能、比热（无混合热）
   3. Arrhenius混合规则计算粘度（粘度对数摩尔分数加权）
   4. Vredeveled混合规则计算热导率（$n=-2$）
7. 向 `equationOfState` 字典增加可选项 `phase` ，允许预测液体密度。
8. 状态方程类实现`dCpdT` `dCvdT`（剩余比热对温度的导数），以供化学反应求解器使用。
   - 目前只有剩余定容比热的导数是严格的。剩余比热差值对温度的导数暂时忽略（太复杂了）

## 修复意外触发SIGFPE（括号备注编译选项）
1. 修复`flowRateInletVelocity`在写出流场时触发SIGFPE。该误触来自于`unitConversion::toUser(const T& t) const`中的除法。（`Clang DP Opt`）
2. 修复`chemistryModel`在计算反应速率触发SIGFPE。该误触来自于`void Foam::Reaction<ThermoType>::C`，可能是分支逻辑导致。（`Clang DP Release`）
3. 全面修复SIGFPE误触。只要有clang且非Debug，就加入`-ffp-exception-behavior`。这不会降低性能，因为主要瓶颈在于内存带宽而非计算速度。

## 修复Bug
1. 修复输运性质为`AndradeTransport`时不能加载化学反应的`dynamicCode`的bug。补上了缺失的函数模板`operator+`和`operator*`的实现。
2. 修复`foamPostProcess`运行一些指定cellZone的后处理函数时崩溃（segfault）的bug。
   - `3caf09d88b95092a1f4a6047cf498d47d5e86a7a` 意外引入这个bug，本意是优化性能。
   - 目前没有找到完美方案修复，暂时回退到旧逻辑。

## 计划添加的
1. 更多状态方程：Patel-Teja、Martin-Hou等
2. 把多孔介质传热模块拓展到单相多组分
3. 改进`porousMediaFluidSolver`的热非平衡模型，使它在大比表面积、高传热系数时更加稳定

## 现存bug(截至2026-04-01):
1. 算例包含拉格朗日场时，`foamRun`崩溃（测试`test/Lagrangian/boundaries`不通过）
