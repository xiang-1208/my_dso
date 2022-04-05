#include "EnergyFunctionalStructs.h"
#include "FullSystem/HessianBlocks.h"

//@ 从 FrameHessian 中提取数据
void EFFrame::takeData()
{
    prior = data->getPrior().head<8>(); 	// 得到先验状态, 主要是光度仿射变换
    delta = data->get_state_minus_stateZero().head<8>(); // 状态与FEJ零状态之间差
    delta_prior =  (data->get_state() - data->getPriorZero()).head<8>(); // 状态与先验之间的差 //? 可先验是0啊?
   

	assert(data->frameID != -1);

	frameID = data->frameID;  // 所有帧的ID序号
}