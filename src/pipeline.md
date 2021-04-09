# pipeline

- 只写单个模板的做法，图形界面交互或者多输入脚本控制都放在外面。

## 流程

- 输入topo结构，用ovm的topology mesh保存
- 计算生成顺序
- 循环
  - deformation（如果是第一步则忽略）
  - 加入（用ovm的geometry mesh保存）
  - 优化（想办法提取关系）
- 输出
