# Todo

写插值

变形前的所有坐标
变形后表面的坐标

用eigen矩阵写

看interpolation的做法
例子里写：
    igl::biharmonic_coordinates(high.V,high.T,S,k,W);
    V、T是V、F，S是一个vector，存控制点的编号，k理论上三维要取3，但是example里说取2问题也不大，W是输出的矩阵

想到一个严肃的问题：如何保障变形的时候仍然无翻转。

整体形变过程：
    控制：表面顶点到目标点
    分段：表面->内部，都用dirichlet能，判断是否反转
    （实际上用的arap）

slim做法：
    确定好b和bc即可，我打算用dirichlet能量
    需要pre_calc
    slim_precompute(V,F,V_0,sData,igl::SLIMData::EXP_CONFORMAL,b,bc,soft_const_p); V是原始，V_0是初始猜测，soft_const_p取1e6这种大数
    slim_solve(sData,1)

arap调用做法：
    需要pre_calc
    igl::arap_precomputation(V,F,V.cols(),b,arap_data);
    igl::arap_solve(bc,arap_data,U);

