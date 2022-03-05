#include "CoarseInitializer.h"

const float minUseGrad_pixsel = 10;

//@ 对于高层(0层以上)选择梯度最大的位置点
template<int pot>
inline int gridMaxSelection(Eigen::Vector3f* grads, bool* map_out, int w, int h, float THFac)
{
	memset(map_out, 0, sizeof(bool)*w*h);

	int numGood = 0;
	for(int y=1;y<h-pot;y+=pot)  // 每隔一个pot遍历
	{
		for(int x=1;x<w-pot;x+=pot)
		{
			int bestXXID = -1; // gradx 最大
			int bestYYID = -1; // grady 最大
			int bestXYID = -1; // gradx-grady 最大
			int bestYXID = -1; // gradx+grady 最大 

			float bestXX=0, bestYY=0, bestXY=0, bestYX=0;

			Eigen::Vector3f* grads0 = grads+x+y*w; // 当前网格的起点
			// 分别找到该网格内上面4个best
			for(int dx=0;dx<pot;dx++)
				for(int dy=0;dy<pot;dy++)
				{
					int idx = dx+dy*w;
					Eigen::Vector3f g=grads0[idx]; // 遍历网格内的每一个像素
					float sqgd = g.tail<2>().squaredNorm(); // 梯度平方和
					float TH = THFac*minUseGrad_pixsel * (0.75f);  //阈值, 为什么都乘0.75 ? downweight

					if(sqgd > TH*TH)
					{
						float agx = fabs((float)g[1]);
						if(agx > bestXX) {bestXX=agx; bestXXID=idx;}

						float agy = fabs((float)g[2]);
						if(agy > bestYY) {bestYY=agy; bestYYID=idx;}

						float gxpy = fabs((float)(g[1]-g[2]));
						if(gxpy > bestXY) {bestXY=gxpy; bestXYID=idx;}

						float gxmy = fabs((float)(g[1]+g[2]));
						if(gxmy > bestYX) {bestYX=gxmy; bestYXID=idx;}
					}
				}

			bool* map0 = map_out+x+y*w; // 选出来的像素为TRUE

			// 选上这些最大的像素
			if(bestXXID>=0)
			{
				if(!map0[bestXXID]) // 没有被选
					numGood++;
				map0[bestXXID] = true;

			}
			if(bestYYID>=0)
			{
				if(!map0[bestYYID])
					numGood++;
				map0[bestYYID] = true;

			}
			if(bestXYID>=0)
			{
				if(!map0[bestXYID])
					numGood++;
				map0[bestXYID] = true;

			}
			if(bestYXID>=0)
			{
				if(!map0[bestYXID])
					numGood++;
				map0[bestYXID] = true;

			}
		}
	}

	return numGood;	
}

//* 同上, 只是把pot作为参数
inline int gridMaxSelection(Eigen::Vector3f* grads, bool* map_out, int w, int h, int pot, float THFac)
{

	memset(map_out, 0, sizeof(bool)*w*h);

	int numGood = 0;
	for(int y=1;y<h-pot;y+=pot)
	{
		for(int x=1;x<w-pot;x+=pot)
		{
			int bestXXID = -1;
			int bestYYID = -1;
			int bestXYID = -1;
			int bestYXID = -1;

			float bestXX=0, bestYY=0, bestXY=0, bestYX=0;

			Eigen::Vector3f* grads0 = grads+x+y*w;
			for(int dx=0;dx<pot;dx++)
				for(int dy=0;dy<pot;dy++)
				{
					int idx = dx+dy*w;
					Eigen::Vector3f g=grads0[idx];
					float sqgd = g.tail<2>().squaredNorm();
					float TH = THFac*minUseGrad_pixsel * (0.75f);

					if(sqgd > TH*TH)
					{
						float agx = fabs((float)g[1]);
						if(agx > bestXX) {bestXX=agx; bestXXID=idx;}

						float agy = fabs((float)g[2]);
						if(agy > bestYY) {bestYY=agy; bestYYID=idx;}

						float gxpy = fabs((float)(g[1]-g[2]));
						if(gxpy > bestXY) {bestXY=gxpy; bestXYID=idx;}

						float gxmy = fabs((float)(g[1]+g[2]));
						if(gxmy > bestYX) {bestYX=gxmy; bestYXID=idx;}
					}
				}

			bool* map0 = map_out+x+y*w;

			if(bestXXID>=0)
			{
				if(!map0[bestXXID])
					numGood++;
				map0[bestXXID] = true;

			}
			if(bestYYID>=0)
			{
				if(!map0[bestYYID])
					numGood++;
				map0[bestYYID] = true;

			}
			if(bestXYID>=0)
			{
				if(!map0[bestXYID])
					numGood++;
				map0[bestXYID] = true;

			}
			if(bestYXID>=0)
			{
				if(!map0[bestYXID])
					numGood++;
				map0[bestYXID] = true;

			}
		}
	}

	return numGood;
}

/********************************
 * @ function:
 * 
 * @ param: 	gards			帧Hessian中当前金字塔层的图像导数
 * @			map				选出的地图点
 * @			desiredDensity	每一金字塔层要的点数(密度)
 * @
 * @ note:		使用递归
 *******************************/
inline int makePixelStatus(Eigen::Vector3f* grads, bool* map, int w, int h, float desiredDensity, int recsLeft=5, float THFac = 1)
{
	if(sparsityFactor < 1) sparsityFactor = 1; // 网格的大小, 在网格内选择最大的

	int numGoodPoints;

	if(sparsityFactor==1) numGoodPoints = gridMaxSelection<1>(grads, map, w, h, THFac);
	else if(sparsityFactor==2) numGoodPoints = gridMaxSelection<2>(grads, map, w, h, THFac);
	else if(sparsityFactor==3) numGoodPoints = gridMaxSelection<3>(grads, map, w, h, THFac);
	else if(sparsityFactor==4) numGoodPoints = gridMaxSelection<4>(grads, map, w, h, THFac);
	else if(sparsityFactor==5) numGoodPoints = gridMaxSelection<5>(grads, map, w, h, THFac);
	else if(sparsityFactor==6) numGoodPoints = gridMaxSelection<6>(grads, map, w, h, THFac);
	else if(sparsityFactor==7) numGoodPoints = gridMaxSelection<7>(grads, map, w, h, THFac);
	else if(sparsityFactor==8) numGoodPoints = gridMaxSelection<8>(grads, map, w, h, THFac);
	else if(sparsityFactor==9) numGoodPoints = gridMaxSelection<9>(grads, map, w, h, THFac);
	else if(sparsityFactor==10) numGoodPoints = gridMaxSelection<10>(grads, map, w, h, THFac);
	else if(sparsityFactor==11) numGoodPoints = gridMaxSelection<11>(grads, map, w, h, THFac);
	else numGoodPoints = gridMaxSelection(grads, map, w, h, sparsityFactor, THFac);

	/*
	 * #points is approximately proportional to sparsityFactor^2.
	 */

	float quotia = numGoodPoints / (float)(desiredDensity);

	int newSparsity = (sparsityFactor * sqrtf(quotia))+0.7f; // 更新网格大小


	if(newSparsity < 1) newSparsity=1;


	float oldTHFac = THFac;
	if(newSparsity==1 && sparsityFactor==1) THFac = 0.5;  // 已经是最小的了, 但是数目还是不够, 就减小阈值

	// 如果满足网格大小变化小且阈值是0.5 || 点数量在20%误差内 || 递归次数已到 , 则返回
	if((abs(newSparsity-sparsityFactor) < 1 && THFac==oldTHFac) ||
			( quotia > 0.8 &&  1.0f / quotia > 0.8) ||
			recsLeft == 0) 
	{

//		printf(" \n");
		//all good
		sparsityFactor = newSparsity;
		return numGoodPoints;
	}
	else // 否则进行递归
	{
//		printf(" -> re-evaluate! \n");
		// re-evaluate.
		sparsityFactor = newSparsity;
		return makePixelStatus(grads, map, w,h, desiredDensity, recsLeft-1, THFac);
	}	
}



CoarseInitializer::CoarseInitializer(int ww, int hh):thisToNext_aff(0,0),thisToNext(SE3())
{
    for(int lvl=0; lvl<pyrLevelsUsed; lvl++)
    {
		points[lvl] = 0;
		numPoints[lvl] = 0;        
    }

    frameID=-1;
}

void CoarseInitializer::setFirst(CalibHessian* HCalib, FrameHessian* newFrameHessian)
{
	//[ ***step 1*** ] 计算图像每层的内参
    makeK(HCalib);
    firstFrame = newFrameHessian;

    PixelSelector sel(w[0],h[0]);

	float* statusMap = new float[w[0]*h[0]];
	bool* statusMapB = new bool[w[0]*h[0]];    

    // 密度权重
    float densities[] = {0.03,0.05,0.15,0.5,1};
    for(int lvl=0; lvl<pyrLevelsUsed; lvl++)
    {
		//[ ***step 2*** ] 针对不同层数选择大梯度像素, 第0层比较复杂1d, 2d, 4d大小block来选择3个层次的像素
        sel.currentPotential = 3;
        int npts;
        if(lvl == 0)	// 第0层提取特征像素
            npts = sel.makeMaps(firstFrame, statusMap,densities[lvl]*w[0]*h[0],1,false,2);
		else	// 其它层则选出goodpoints
			npts = makePixelStatus(firstFrame->dIp[lvl], statusMapB, w[lvl], h[lvl], densities[lvl]*w[0]*h[0]);

		// 如果点非空, 则释放空间, 创建新的
		if(points[lvl] != 0) delete[] points[lvl];
		points[lvl] = new Pnt[npts];

		int wl=w[lvl], hl=h[lvl]; 
		Pnt* pl = points[lvl];
		int nl = 0;

		// 要留出pattern的空间, 2 border
		//[ ***step 3*** ] 在选出的像素中, 添加点信息
		for(int y=patternPadding+1;y<hl-patternPadding-2;y++)
		for(int x=patternPadding+1;x<wl-patternPadding-2;x++)
		{
			if((lvl!=0 && statusMapB[x+y*wl]) || (lvl==0 && statusMap[x+y*wl] != 0))
			{
				//assert(patternNum==9);
				pl[nl].u = x+0.1;   //? 加0.1干啥
				pl[nl].v = y+0.1;
				pl[nl].idepth = 1;
				pl[nl].iR = 1;
				pl[nl].isGood=true;
				pl[nl].energy.setZero();
				pl[nl].lastHessian=0;
				pl[nl].lastHessian_new=0;
				pl[nl].my_type= (lvl!=0) ? 1 : statusMap[x+y*wl];

				Eigen::Vector3f* cpt = firstFrame->dIp[lvl] +x + y*w[lvl];	// 该像素梯度
				float sumGrad2=0;
				// 计算pattern内像素梯度和
				for(int idx=0;idx<patternNum;idx++)
				{
					int dx = patternP[idx][0]; // pattern 的偏移
					int dy = patternP[idx][1];
					float absgrad = cpt[dx + dy*w[lvl]].tail<2>().squaredNorm();
					sumGrad2 += absgrad;
				}

				// float gth = setting_outlierTH * (sqrtf(sumGrad2)+setting_outlierTHSumComponent);
				// pl[nl].outlierTH = patternNum*gth*gth;
				//! 外点的阈值与pattern的大小有关, 一个像素是12*12
				//? 这个阈值怎么确定的...
				pl[nl].outlierTH = patternNum*setting_outlierTH;		

				nl++;
				assert(nl <= npts);
			}
		}	

		numPoints[lvl]=nl; // 点的数目,  去掉了一些边界上的点	
    }
	delete[] statusMap;
	delete[] statusMapB;	

	//[ ***step 4*** ] 计算点的最近邻和父点
	makeNN();

	// 参数初始化
	thisToNext=SE3();
	snapped = false;
	frameID = snappedAt = 0;


	for(int i=0;i<pyrLevelsUsed;i++)
		dGrads[i].setZero();	
}	

void CoarseInitializer::makeK(CalibHessian* HCalib)
{
	w[0] = wG[0];
	h[0] = hG[0];

	fx[0] = HCalib->fxl();
	fy[0] = HCalib->fyl();
	cx[0] = HCalib->cxl();
	cy[0] = HCalib->cyl();    

    for (int level = 1; level < pyrLevelsUsed; ++ level)
    {
		w[level] = w[0] >> level;
		h[level] = h[0] >> level;
		fx[level] = fx[level-1] * 0.5;
		fy[level] = fy[level-1] * 0.5;
		cx[level] = (cx[0] + 0.5) / ((int)1<<level) - 0.5;
		cy[level] = (cy[0] + 0.5) / ((int)1<<level) - 0.5;        
    }

	for (int level = 0; level < pyrLevelsUsed; ++ level)
	{
		K[level]  << fx[level], 0.0, cx[level], 0.0, fy[level], cy[level], 0.0, 0.0, 1.0;
		Ki[level] = K[level].inverse();
		fxi[level] = Ki[level](0,0);
		fyi[level] = Ki[level](1,1);
		cxi[level] = Ki[level](0,2);
		cyi[level] = Ki[level](1,2);
	}
}

//@ 生成每一层点的KDTree, 并用其找到邻近点集和父点 
void CoarseInitializer::makeNN()
{
	const float NNDistFactor=0.05;
	// 第一个参数为distance, 第二个是datasetadaptor, 第三个是维数
	typedef nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<float, FLANNPointcloud> ,
			FLANNPointcloud,2> KDTree;

	// build indices
	FLANNPointcloud pcs[PYR_LEVELS]; // 每层建立一个点云	
	KDTree* indexes[PYR_LEVELS]; // 点云建立KDtree
	//* 每层建立一个KDTree索引二维点云
	for(int i=0;i<pyrLevelsUsed;i++)
	{
		pcs[i] = FLANNPointcloud(numPoints[i], points[i]); // 二维点点云
		// 参数: 维度, 点数据, 叶节点中最大的点数(越大build快, query慢)
		indexes[i] = new KDTree(2, pcs[i], nanoflann::KDTreeSingleIndexAdaptorParams(5) );
		indexes[i]->buildIndex();
	}

	const int nn=10;

	// find NN & parents
	for(int lvl=0;lvl<pyrLevelsUsed;lvl++)
	{
		Pnt* pts = points[lvl];
		int npts = numPoints[lvl];

		int ret_index[nn];  // 搜索到的临近点
		float ret_dist[nn]; // 搜索到点的距离
		// 搜索结果, 最近的nn个和1个
		nanoflann::KNNResultSet<float, int, int> resultSet(nn);
		nanoflann::KNNResultSet<float, int, int> resultSet1(1);

		for(int i=0;i<npts;i++)
		{
			//resultSet.init(pts[i].neighbours, pts[i].neighboursDist );
			resultSet.init(ret_index, ret_dist);
			Vec2f pt = Vec2f(pts[i].u,pts[i].v); // 当前点
			// 使用建立的KDtree, 来查询最近邻
			indexes[lvl]->findNeighbors(resultSet, (float*)&pt, nanoflann::SearchParams());
			int myidx=0;
			float sumDF = 0;
			
			//* 给每个点的neighbours赋值
			for(int k=0;k<nn;k++)
			{
				pts[i].neighbours[myidx]=ret_index[k]; // 最近的索引
				float df = expf(-ret_dist[k]*NNDistFactor); // 距离使用指数形式
				sumDF += df; // 距离和
				pts[i].neighboursDist[myidx]=df;
				assert(ret_index[k]>=0 && ret_index[k] < npts);
				myidx++;
			}

			// 对距离进行归10化,,,,,
			for(int k=0;k<nn;k++)
				pts[i].neighboursDist[k] *= 10/sumDF;	


			//* 高一层的图像中找到该点的父节点
			if(lvl < pyrLevelsUsed-1 )	
			{
				resultSet1.init(ret_index, ret_dist);
				pt = pt*0.5f-Vec2f(0.25f,0.25f); // 换算到高一层
				indexes[lvl+1]->findNeighbors(resultSet1, (float*)&pt, nanoflann::SearchParams());
				pts[i].parent = ret_index[0]; // 父节点
				pts[i].parentDist = expf(-ret_dist[0]*NNDistFactor); // 到父节点的距离(在高层中)

				assert(ret_index[0]>=0 && ret_index[0] < numPoints[lvl+1]);
			}
			else  // 最高层没有父节点
			{
				pts[i].parent = -1;
				pts[i].parentDist = -1;
			}
		}
	}	
}

bool CoarseInitializer::trackFrame(FrameHessian* newFrameHessian, std::vector<IOWrap::Output3DWrapper*> &wraps)
{
	newFrame = newFrameHessian;
	std::cout<<"!!!!!!"<<std::endl;
	//[ ***step 1*** ] 先显示新来的帧
	// 新的一帧, 在跟踪之前显示的
	for(IOWrap::Output3DWrapper* ow : wraps)
        ow->pushLiveFrame(newFrameHessian);

	int maxIterations[] = {5,5,10,30,50};

	//????? 调参
	alphaK = 2.5*2.5;//*freeDebugParam1*freeDebugParam1;
	alphaW = 150*150;//*freeDebugParam2*freeDebugParam2;
	regWeight = 0.8;//*freeDebugParam4;
	couplingWeight = 1;//*freeDebugParam5;	

	//[ ***step 2*** ] 初始化每个点逆深度为1, 初始化光度参数, 位姿SE3
	if(!snapped) 
	{
		// 初始化
		thisToNext.translation().setZero();
		for(int lvl=0;lvl<pyrLevelsUsed;lvl++)
		{
			int npts = numPoints[lvl];
			Pnt* ptsl = points[lvl];
			for (int i=0;i<npts;i++)
			{
				ptsl[i].iR = 1;
				ptsl[i].idepth_new = 1;
				ptsl[i].lastHessian = 0;				
			}
		}
	}

	SE3 refToNew_current = thisToNext;
	AffLight refToNew_aff_current = thisToNext_aff;

	// 如果都有仿射系数, 则估计一个初值
	if(firstFrame->ab_exposure>0 && newFrame->ab_exposure>0)
		refToNew_aff_current = AffLight(logf(newFrame->ab_exposure /  firstFrame->ab_exposure),0); // coarse approximation.

	Vec3f latestRes = Vec3f::Zero();
	// 从顶层开始估计
	for(int lvl=pyrLevelsUsed-1; lvl>=0; lvl--)
	{
		//[ ***step 3*** ] 使用计算过的上一层来初始化下一层
		// 顶层未初始化到, reset来完成
		if(lvl<pyrLevelsUsed-1)
			propagateDown(lvl+1);		

		Mat88f H,Hsc; Vec8f b,bsc;
		resetPoints(lvl); //这里对顶层进行初始化！
		//[ ***step 4*** ] 迭代之前计算能量, Hessian等
		Vec3f resOld = calcResAndGS(lvl, H, b, Hsc, bsc, refToNew_current, refToNew_aff_current, false);

	}

	return true;
}

//@ 使用上层信息来初始化下层
//@ param: 当前的金字塔层+1
//@ note: 没法初始化顶层值 
void CoarseInitializer::propagateDown(int srcLvl)
{
	assert(srcLvl>0);

	int nptst= numPoints[srcLvl-1]; // 当前层的点数目
	Pnt* ptss = points[srcLvl];  // 当前层+1, 上一层的点集
	Pnt* ptst = points[srcLvl-1]; // 当前层点集

	for(int i=0;i<nptst;i++)
	{
		Pnt* point = ptst+i;  // 遍历当前层的点
		Pnt* parent = ptss+point->parent;  // 找到当前点的parrent

		if(!parent->isGood || parent->lastHessian < 0.1) continue;
		if(!point->isGood)
		{
			// 当前点不好, 则把父点的值直接给它, 并且置位good
			point->iR = point->idepth = point->idepth_new = parent->iR;
			point->isGood=true;
			point->lastHessian=0;
		}
		else
		{
			// 通过hessian给point和parent加权求得新的iR
			// iR可以看做是深度的值, 使用的高斯归一化积, Hessian是信息矩阵	
			//??????
			float newiR = (point->iR*point->lastHessian*2 + parent->iR*parent->lastHessian) / (point->lastHessian*2+parent->lastHessian);
			point->iR = point->idepth = point->idepth_new = newiR;		
		}
	}
	//? 为什么在这里又更新了iR, 没有更新 idepth 
	// 感觉更多的是考虑附近点的平滑效果
	optReg(srcLvl-1); // 当前层
}

//* 使用最近点来更新每个点的iR, smooth的感觉
void CoarseInitializer::optReg(int lvl)
{
	int npts = numPoints[lvl];
	Pnt* ptsl = points[lvl];

	//* 位移不足够则设置iR是1
	if(!snapped)
	{
		for(int i=0;i<npts;i++)
			ptsl[i].iR = 1;
		return;
	}

	for(int i=0;i<npts;i++)
	{
		Pnt* point = ptsl+i;
		if(!point->isGood) continue;

		float idnn[10];
		int nnn=0;
		// 获得当前点周围最近10个点, 质量好的点的iR
		for(int j=0;j<10;j++)
		{
			if(point->neighbours[j] == -1) continue;
			Pnt* other = ptsl+point->neighbours[j];
			if(!other->isGood) continue;
			idnn[nnn] = other->iR;
			nnn++;
		}

		// 与最近点中位数进行加权获得新的iR
		if(nnn > 2)
		{
			std::nth_element(idnn,idnn+nnn/2,idnn+nnn); // 获得中位数
			point->iR = (1-regWeight)*point->idepth + regWeight*idnn[nnn/2];
		}
	}
}

//@ 重置点的energy, idepth_new参数
void CoarseInitializer::resetPoints(int lvl)
{
	Pnt* pts = points[lvl];
	int npts = numPoints[lvl];
	for (int i=0;i<npts;i++)
	{
		//重置
		pts[i].energy.setZero();
		pts[i].idepth_new = pts[i].idepth;

		// 如果是最顶层, 则使用周围点平均值来重置
		if(lvl==pyrLevelsUsed-1 && !pts[i].isGood)
		{
			float snd=0,sn=0;
			for(int n=0;n<10;n++)
			{
				if(pts[i].neighbours[n] == -1 || !pts[pts[i].neighbours[n]].isGood) continue;
				snd += pts[pts[i].neighbours[n]].iR;
				sn += 1;
			}

			if(sn > 0)
			{
				pts[i].isGood=true;
				pts[i].iR = pts[i].idepth = pts[i].idepth_new = snd / sn;
			}
		}
	}
}

//* 计算能量函数和Hessian矩阵, 以及舒尔补, sc代表Schur
// calculates residual, Hessian and Hessian-block neede for re-substituting depth.
Vec3f CoarseInitializer::calcResAndGS(
		int lvl, Mat88f &H_out, Vec8f &b_out,
		Mat88f &H_out_sc, Vec8f &b_out_sc,
		const SE3 &refToNew, AffLight refToNew_aff,
		bool plot)
{
	int wl = w[lvl], hl = h[lvl];
	// 当前层图像及梯度
	Eigen::Vector3f* colorRef = firstFrame->dIp[lvl];  
	Eigen::Vector3f* colorNew = newFrame->dIp[lvl];

	//! 旋转矩阵R * 内参矩阵K_inv
	Mat33f RKi = (refToNew.rotationMatrix() * Ki[lvl]).cast<float>();	
	Vec3f t = refToNew.translation().cast<float>(); // 平移
	Eigen::Vector2f r2new_aff = Eigen::Vector2f(exp(refToNew_aff.a), refToNew_aff.b); // 光度参数

	// 该层的相机参数
	float fxl = fx[lvl];
	float fyl = fy[lvl];
	float cxl = cx[lvl];
	float cyl = cy[lvl];


	Accumulator11 E;  // 1*1 的累加器
	acc9.initialize(); // 初始值, 分配空间
	E.initialize();


	int npts = numPoints[lvl];
	Pnt* ptsl = points[lvl];
	for(int i=0;i<npts;i++)
	{

		Pnt* point = ptsl+i;

		point->maxstep = 1e10;
		if(!point->isGood)  // 点不好
		{
			E.updateSingle((float)(point->energy[0])); // 累加
			point->energy_new = point->energy;
			point->isGood_new = false;
			continue;
		}

        VecNRf dp0;  // 8*1矩阵, 每个点附近的残差个数为8个
        VecNRf dp1;
        VecNRf dp2;
        VecNRf dp3;
        VecNRf dp4;
        VecNRf dp5;
        VecNRf dp6;
        VecNRf dp7;
        VecNRf dd;
        VecNRf r;
		JbBuffer_new[i].setZero();  // 10*1 向量

		// sum over all residuals.
		bool isGood = true;
		float energy=0;
		for(int idx=0;idx<patternNum;idx++)
		{
			// pattern的坐标偏移
			int dx = patternP[idx][0];
			int dy = patternP[idx][1];

			//! Pj' = R*(X/Z, Y/Z, 1) + t/Z, 变换到新的点, 深度仍然使用Host帧的!
			Vec3f pt = RKi * Vec3f(point->u+dx, point->v+dy, 1) + t*point->idepth_new; 
			// 归一化坐标 Pj
			float u = pt[0] / pt[2];
			float v = pt[1] / pt[2];
			// 像素坐标pj
			float Ku = fxl * u + cxl;
			float Kv = fyl * v + cyl;
			// dpi/pz' 
			float new_idepth = point->idepth_new/pt[2]; // 新一帧上的逆深度

			// 落在边缘附近，深度小于0, 则不好
			if(!(Ku > 1 && Kv > 1 && Ku < wl-2 && Kv < hl-2 && new_idepth > 0))
			{
				isGood = false;
				break;
			}
			// 插值得到新图像中的 patch 像素值，(输入3维，输出3维像素值 + x方向梯度 + y方向梯度)
			Vec3f hitColor = getInterpolatedElement33(colorNew, Ku, Kv, wl);
			//Vec3f hitColor = getInterpolatedElement33BiCub(colorNew, Ku, Kv, wl);

			// 参考帧上的 patch 上的像素值, 输出一维像素值
			//float rlR = colorRef[point->u+dx + (point->v+dy) * wl][0];
			float rlR = getInterpolatedElement31(colorRef, point->u+dx, point->v+dy, wl);

			// 像素值有穷, good
			if(!std::isfinite(rlR) || !std::isfinite((float)hitColor[0]))
			{
				isGood = false;
				break;
			}

			// 残差
			float residual = hitColor[0] - r2new_aff[0] * rlR - r2new_aff[1];
			// Huber权重
			float hw = fabs(residual) < setting_huberTH ? 1 : setting_huberTH / fabs(residual); 
			// huberweight * (2-huberweight) = Objective Function
			// robust 权重和函数之间的关系
			energy += hw *residual*residual*(2-hw);

			// Pj 对 逆深度 di 求导   
			//! 1/Pz * (tx - u*tz), u = px/pz 
			float dxdd = (t[0]-t[2]*u)/pt[2];   
			//! 1/Pz * (ty - v*tz), u = py/pz
			float dydd = (t[1]-t[2]*v)/pt[2];

			// 参考https://www.cnblogs.com/JingeTU/p/8203606.html
			if(hw < 1) hw = sqrtf(hw); //?? 为啥开根号, 答: 鲁棒核函数等价于加权最小二乘
			//! dxfx, dyfy
			float dxInterp = hw*hitColor[1]*fxl;
			float dyInterp = hw*hitColor[2]*fyl;
			//* 残差对 j(新状态) 位姿求导, 
			dp0[idx] = new_idepth*dxInterp; //! dpi/pz' * dxfx
			dp1[idx] = new_idepth*dyInterp; //! dpi/pz' * dyfy
			dp2[idx] = -new_idepth*(u*dxInterp + v*dyInterp); //! -dpi/pz' * (px'/pz'*dxfx + py'/pz'*dyfy)
			dp3[idx] = -u*v*dxInterp - (1+v*v)*dyInterp; //! - px'py'/pz'^2*dxfy - (1+py'^2/pz'^2)*dyfy
			dp4[idx] = (1+u*u)*dxInterp + u*v*dyInterp; //! (1+px'^2/pz'^2)*dxfx + px'py'/pz'^2*dxfy
			dp5[idx] = -v*dxInterp + u*dyInterp; //! -py'/pz'*dxfx + px'/pz'*dyfy
			//* 残差对光度参数求导
			dp6[idx] = - hw*r2new_aff[0] * rlR; //! exp(aj-ai)*I(pi)
			dp7[idx] = - hw*1;	//! 对 b 导
			//* 残差对 i(旧状态) 逆深度求导
			dd[idx] = dxInterp * dxdd  + dyInterp * dydd; 	//! dxfx * 1/Pz * (tx - u*tz) +　dyfy * 1/Pz * (tx - u*tz)			
			r[idx] = hw*residual; //! 残差 res

			
		}
	}
}