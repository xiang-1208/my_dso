#pragma once
#include <Eigen/Core>
#include "FullSystem/FrameShell.h"
#include "util/globalutil.h"


using namespace dso;

class EFFrame;
struct FrameHessian;
struct CalibHessian;
struct PointHessian;

//? 这是干什么用的? 是为了求解时候的数值稳定? 
#define SCALE_IDEPTH 1.0f			//!< 逆深度的比例系数  // scales internal value to idepth.
#define SCALE_XI_ROT 1.0f			//!< 旋转量(so3)的比例系数
#define SCALE_XI_TRANS 0.5f			//!< 平移量的比例系数, 尺度?
#define SCALE_F 50.0f   			//!< 相机焦距的比例系数
#define SCALE_C 50.0f				//!< 相机光心偏移的比例系数
#define SCALE_W 1.0f				//!< 不知道...
#define SCALE_A 10.0f				//!< 光度仿射系数a的比例系数
#define SCALE_B 1000.0f				//!< 光度仿射系数b的比例系数

//* 其中带0的是FEJ用的初始状态, 不带0的是更新的状态
struct FrameFramePrecalc
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	FrameHessian* host;	// defines row
	FrameHessian* target;	// defines column

	Mat33f PRE_RTll; // host 到 target 之间优化后旋转矩阵 R
	Mat33f PRE_KRKiTll; // k*R*k_inv
	Mat33f PRE_RKiTll;  // R*k_inv
	Mat33f PRE_RTll_0; // host 到 target之间初始的旋转矩阵, 优化更新前

	Vec3f PRE_tTll; //  host 到 target之间优化后的平移 t
	Vec3f PRE_KtTll; // K*t
	Vec3f PRE_tTll_0; //  host 到 target之间初始的平移, 优化更新前

	float distanceLL; // 两帧间距离

	Vec2f PRE_aff_mode; // 能量函数对仿射系数处理后的, 总系数
	float PRE_b0_mode; // host的光度仿射系数b

	void set(FrameHessian* host, FrameHessian* target, CalibHessian* HCalib);	
};

struct CalibHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static int instanceCounter;

    float Binv[256];
    float B[256];

    VecCf value_scaledf;
	VecC value_minus_value_zero;	//!< 减去线性化点

    inline ~CalibHessian() {instanceCounter--;}
	inline CalibHessian()
	{
		VecC initial_value = VecC::Zero();
		initial_value[0] = fxG[0];
		initial_value[1] = fyG[0];
		initial_value[2] = cxG[0];
		initial_value[3] = cyG[0];

		// setValueScaled(initial_value);
		// value_zero = value;
		value_minus_value_zero.setZero();

		instanceCounter++;
		for(int i=0;i<256;i++)
			Binv[i] = B[i] = i;		// set gamma function to identity
	};

    EIGEN_STRONG_INLINE float getBGradOnly(float color)
	{
		int c = color+0.5f;
		if(c<5) c=5;
		if(c>250) c=250;
		return B[c+1]-B[c];
	}

    inline float& fxl() {return value_scaledf[0];}
    inline float& fyl() {return value_scaledf[1];}
    inline float& cxl() {return value_scaledf[2];}
    inline float& cyl() {return value_scaledf[3];}
};

struct FrameHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EFFrame* efFrame;		//!< 帧的能量函数

    FrameShell* shell;
    float ab_exposure;

	std::vector<PointHessian*> pointHessians;				//!< contains all ACTIVE points.
	std::vector<PointHessian*> pointHessiansMarginalized;	//!< contains all MARGINALIZED points (= fully marginalized, usually because point went OOB.)
	std::vector<PointHessian*> pointHessiansOut;			//!< contains all OUTLIER points (= discarded.).
	//std::vector<ImmaturePoint*> immaturePoints;				//!< contains all OUTLIER points (= discarded.).

    Eigen::Vector3f* dI;				//!< 图像导数 // trace, fine tracking. Used for direction select (not for gradient histograms etc.)
	Eigen::Vector3f* dIp[PYR_LEVELS];	 // coarse tracking / coarse initializer. NAN in [0] only.
    float* absSquaredGrad[PYR_LEVELS];

	SE3 PRE_worldToCam;			//!< 预计算的, 位姿状态增量更新到位姿上
	SE3 PRE_camToWorld;
	std::vector<FrameFramePrecalc,Eigen::aligned_allocator<FrameFramePrecalc>> targetPrecalc; //!< 对于其它帧的预运算值

	int frameID;					//!< 所有关键帧的序号(FrameShell)
	int idx;						//!< 激活关键帧的序号(FrameHessian)

	SE3 worldToCam_evalPT;		//!< 在估计的相机位姿
	//* 这三个是与线性化点的增量, 而光度参数不是增量, state就是值
	Vec10 state_zero;   		//!< 固定的线性化点的状态增量, 为了计算进行缩放
	Vec10 state_scaled;			//!< 乘上比例系数的状态增量, 这个是真正求的值!!!
	Vec10 state;				//!< 计算的状态增量	

	inline FrameHessian()
	{
		frameID = -1;
	};

	EIGEN_STRONG_INLINE const SE3 &get_worldToCam_evalPT() const {return worldToCam_evalPT;}
    EIGEN_STRONG_INLINE const Vec10 &get_state_zero() const {return state_zero;}
    EIGEN_STRONG_INLINE const Vec10 &get_state() const {return state;}
	EIGEN_STRONG_INLINE const Vec10 &get_state_scaled() const {return state_scaled;}
	EIGEN_STRONG_INLINE const Vec10 get_state_minus_stateZero() const {return get_state()-get_state_zero();}

	//* 获得先验信息矩阵
	inline Vec10 getPrior()
	{
		Vec10 p =  Vec10::Zero();

		if(frameID==0)  //* 第一帧就用初始值做先验
		{
			p.head<3>() = Vec3::Constant(setting_initialTransPrior);
			p.segment<3>(3) = Vec3::Constant(setting_initialRotPrior);


			p[6] = setting_initialAffAPrior; // 1e14
			p[7] = setting_initialAffBPrior; // 1e14						
		}
		else
		{
			if(setting_affineOptModeA < 0) //* 小于零是固定的不优化
				p[6] = setting_initialAffAPrior;
			else
				p[6] = setting_affineOptModeA;   // 1e12

			if(setting_affineOptModeB < 0)
				p[7] = setting_initialAffBPrior;
			else
				p[7] = setting_affineOptModeB;  // 1e8			
		}
		//? 8,9是干嘛的呢???  没用....
		p[8] = setting_initialAffAPrior;
		p[9] = setting_initialAffBPrior;
		return p;
	}

	inline AffLight aff_g2l() const {return AffLight(get_state_scaled()[6], get_state_scaled()[7]);} //* 返回光度仿射系数
	inline AffLight aff_g2l_0()	const {return AffLight(get_state_zero()[6]*SCALE_A, get_state_zero()[7]*SCALE_B);} //* 返回线性化点处的仿射系数增量

    void makeImages(float* color, CalibHessian* HCalib);

	inline Vec10 getPriorZero()
	{
		return Vec10::Zero();
	}
};

//* 点Hessian
// hessian component associated with one point.
struct PointHessian
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	float idepth_zero;					//!< 缩放了scale倍的固定线性化点逆深度
	float idepth;						//!< 缩放scale倍的逆深度	

	PointHessian(const ImmaturePoint* const rawPoint, CalibHessian* Hcalib);
};

