#include "PixelSelector2.h"


int computeHistQuantil(int* hist, float below)
{
	int th = hist[0]*below+0.5f;
	for(int i=0;i<90;i++)
	{
		th -= hist[i+1];
		if(th<0) return i;
	}
	return 90;                                 
}

PixelSelector::PixelSelector(int w, int h)
{
    randomPattern = new unsigned char[w*h];
    std::srand(3141592);
    for(int i=0;i<w*h;i++) randomPattern[i] = rand() & 0xFF;

    currentPotential=3;

    gradHistFrame=0;
    gradHist = new int[100*(1+w/32)*(1+h/32)];
    ths = new float[(w/32)*(h/32)+100];
    thsSmoothed = new float[(w/32)*(h/32)+100];
}


/********************************
 * @ function:
 * 
 * @ param: 	fh				帧Hessian数据结构
 * @			map_out			选出的地图点
 * @			density		 	每一金字塔层要的点数(密度)
 * @			recursionsLeft	最大递归次数
 * @			plot			画图
 * @			thFactor		阈值因子
 * @
 * @ note:		使用递归
 *******************************/
int PixelSelector::makeMaps(
		const FrameHessian* const fh,
		float* map_out, float density, int recursionsLeft, bool plot, float thFactor)
{
	float numHave=0;
	float numWant=density;
	float quotia;
	int idealPotential = currentPotential;  

	//[ ***step 1*** ] 没有计算直方图, 以及选点的阈值, 则调用函数生成block阈值
    if(fh != gradHistFrame) makeHists(fh);  

    //[ ***step 2*** ] 在当前帧上选择符合条件的像素
    Eigen::Vector3i n = this->select(fh, map_out,currentPotential, thFactor);

	numHave = n[0]+n[1]+n[2]; // 选择得到的点
	quotia = numWant / numHave;  // 得到的 与 想要的 比例

	//[ ***step 3*** ] 计算新的采像素点的, 范围大小, 相当于动态网格了, pot越小取得点越多
	// by default we want to over-sample by 40% just to be sure.
	float K = numHave * (currentPotential+1) * (currentPotential+1);	// 相当于覆盖的面积, 每一个像素对应一个pot*pot
	idealPotential = sqrtf(K/numWant)-1;	// round down.
	if(idealPotential<1) idealPotential=1;

	//[ ***step 4*** ] 想要的数目和已经得到的数目, 大于或小于0.25都会重新采样一次
	if( recursionsLeft>0 && quotia > 1.25 && currentPotential>1)
	{
		if(idealPotential>=currentPotential)	// idealPotential应该小
			idealPotential = currentPotential-1;	// 减小,多采点
		currentPotential = idealPotential;		
		return makeMaps(fh,map_out, density, recursionsLeft-1, plot,thFactor); //递归
	}
	else if(recursionsLeft>0 && quotia < 0.25)
	{
		if(idealPotential<=currentPotential)	// idealPotential应该大
			idealPotential = currentPotential+1;	// 增大, 少采点
		currentPotential = idealPotential;		
		return makeMaps(fh,map_out, density, recursionsLeft-1, plot,thFactor); //递归		
	}

	//[ ***step 5*** ] 现在提取的还是多, 随机删除一些点
	int numHaveSub = numHave;
	if(quotia < 0.95)
	{
		int wh=wG[0]*hG[0];
		int rn=0;
		unsigned char charTH = 255*quotia;
		for(int i=0;i<wh;i++)
		{
			if(map_out[i] != 0)
			{
				if(randomPattern[rn] > charTH )
				{
					map_out[i]=0;
					numHaveSub--;
				}
				rn++;
			}
		}
	}

	currentPotential = idealPotential;


	// 画出选择结果
	if(plot)
	{
		int w = wG[0];
		int h = hG[0];


		MinimalImageB3 img(w,h);

		for(int i=0;i<w*h;i++)
		{
			float c = fh->dI[i][0]*0.7; // 像素值
			if(c>255) c=255;
			img.at(i) = Vec3b(c,c,c);
		}
		IOWrap::displayImage("Selector Image", &img);

		// 安照不同层数的像素, 画上不同颜色
		for(int y=0; y<h;y++)
			for(int x=0;x<w;x++)
			{
				int i=x+y*w;
				if(map_out[i] == 1)
					img.setPixelCirc(x,y,Vec3b(0,255,0));
				else if(map_out[i] == 2)
					img.setPixelCirc(x,y,Vec3b(255,0,0));
				else if(map_out[i] == 4)
					img.setPixelCirc(x,y,Vec3b(0,0,255));
			}
		IOWrap::displayImage("Selector Pixels", &img);
	}

	return numHaveSub;	
}

// 计算32*32局部区域块的梯度中值
void PixelSelector::makeHists(const FrameHessian* const fh)
{
    // 要选取特征点的一帧结构
    gradHistFrame = fh;
    // 提前计算的图像梯度（dx*dx + dy*dy）
    float * mapmax0 = fh->absSquaredGrad[0];

	int w = wG[0];
	int h = hG[0];

	int w32 = w/32;
	int h32 = h/32;
	thsStep = w32;


    // 统计32 X 32每个区域块中的梯度直方图, 只要某点的梯度大于48, 就将该点的梯度设置为48
    // 其中hist[0]统计区域块中总的数目
    // hist[1]至hist[48]分别统计梯度为(0-47)的个数
    // hist[49]统计梯度大于等于48的个数
	for(int y=0;y<h32;y++)
		for(int x=0;x<w32;x++)
        {
            float* map0 = mapmax0+32*x+32*y*w;
            int* hist0 = gradHist;
            memset(hist0,0,sizeof(int)*50);

            for(int j=0;j<32;j++) for(int i=0;i<32;i++)
            {
				int it = i+32*x;
				int jt = j+32*y;
                if(it>w-2 || jt>h-2 || it<1 || jt<1) continue;
                int g = sqrtf(map0[i+j*w]);
                if(g>48) g=48;
				hist0[g+1]++;
				hist0[0]++;
            }

            // 根据设定的阈值(setting_minGradHistCut)来设定区域阈值
            // 比如默认设定setting_minGradHistCut = 0.5f 就是选取中值阈值
            ths[x+y*w32] = computeHistQuantil(hist0,setting_minGradHistCut) + setting_minGradHistAdd;           
        }

    // 这一部分就是对32*32个计算的区域阈值进行3邻域的均值平滑
	for(int y=0;y<h32;y++)
		for(int x=0;x<w32;x++)
        {
            float sum=0,num=0;
            if(x>0)
            {
                if(y>0) 	{num++; 	sum+=ths[x-1+(y-1)*w32];}
                if(y<h32-1) {num++; 	sum+=ths[x-1+(y+1)*w32];}
            }

			if(x<w32-1)
			{
				if(y>0) 	{num++; 	sum+=ths[x+1+(y-1)*w32];}
				if(y<h32-1) {num++; 	sum+=ths[x+1+(y+1)*w32];}
				num++; sum+=ths[x+1+(y)*w32];
            }

			if(y>0) 	{num++; 	sum+=ths[x+(y-1)*w32];}
			if(y<h32-1) {num++; 	sum+=ths[x+(y+1)*w32];}
			num++; sum+=ths[x+y*w32];

			thsSmoothed[x+y*w32] = (sum/num) * (sum/num);
        }    
}

Eigen::Vector3i PixelSelector::select(const FrameHessian* const fh,
		float* map_out, int pot, float thFactor)
{ 
    Eigen::Vector3f const * const map0 = fh->dI;

	float * mapmax0 = fh->absSquaredGrad[0];
	float * mapmax1 = fh->absSquaredGrad[1];
	float * mapmax2 = fh->absSquaredGrad[2];

	int w = wG[0];
	int w1 = wG[1];
	int w2 = wG[2];
	int h = hG[0];   

	const Vec2f directions[16] = {
	         Vec2f(0,    1.0000),
	         Vec2f(0.3827,    0.9239),
	         Vec2f(0.1951,    0.9808),
	         Vec2f(0.9239,    0.3827),
	         Vec2f(0.7071,    0.7071),
	         Vec2f(0.3827,   -0.9239),
	         Vec2f(0.8315,    0.5556),
	         Vec2f(0.8315,   -0.5556),
	         Vec2f(0.5556,   -0.8315),
	         Vec2f(0.9808,    0.1951),
	         Vec2f(0.9239,   -0.3827),
	         Vec2f(0.7071,   -0.7071),
	         Vec2f(0.5556,    0.8315),
	         Vec2f(0.9808,   -0.1951),
	         Vec2f(1.0000,    0.0000),
	         Vec2f(0.1951,   -0.9808)}; 

    memset(map_out,0,w*h*sizeof(PixelSelectorStatus));

	float dw1 = setting_gradDownweightPerLevel;
	float dw2 = dw1*dw1; 

    int n3=0, n2=0, n4=0;  
    // 4d * 4d 区域内选取
    for(int y4=0;y4<h;y4+=(4*pot)) for(int x4=0;x4<w;x4+=(4*pot))
    {
        int my3 = std::min((4*pot), h-y4);
        int mx3 = std::min((4*pot), w-x4);
        int bestIdx4=-1; float bestVal4=0;
        Vec2f dir4 = directions[randomPattern[n2] & 0xF];
        for(int y3=0;y3<my3;y3+=(2*pot)) for(int x3=0;x3<mx3;x3+=(2*pot))
        {
            int x34 = x3+x4;
            int y34 = y3+y4;
            int my2 = std::min((2*pot), h-y34);
			int mx2 = std::min((2*pot), w-x34);
            int bestIdx3=-1; float bestVal3=0;
            Vec2f dir3 = directions[randomPattern[n2] & 0xF];
            for(int y2=0;y2<my2;y2+=pot) for(int x2=0;x2<mx2;x2+=pot)
            {
				int x234 = x2+x34;
				int y234 = y2+y34;
				int my1 = std::min(pot, h-y234);
				int mx1 = std::min(pot, w-x234);
				int bestIdx2=-1; float bestVal2=0;
				Vec2f dir2 = directions[randomPattern[n2] & 0xF];  
                for(int y1=0;y1<my1;y1+=1) for(int x1=0;x1<mx1;x1+=1)
                {
					assert(x1+x234 < w);
					assert(y1+y234 < h); 
                    int idx = x1+x234 + w*(y1+y234);   
                    int xf = x1+x234;
					int yf = y1+y234;    

                    if(xf<4 || xf>=w-5 || yf<4 || yf>h-4) continue;  

					float pixelTH0 = thsSmoothed[(xf>>5) + (yf>>5) * thsStep];
					float pixelTH1 = pixelTH0*dw1;
					float pixelTH2 = pixelTH1*dw2;    

					float ag0 = mapmax0[idx];                    
					if(ag0 > pixelTH0*thFactor)
					{
                        //获取向量尾部的n个元素：vector.tail(n);
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir2)));
						if(!setting_selectDirectionDistribution) dirNorm = ag0;

						if(dirNorm > bestVal2)
						{ 
                            bestVal2 = dirNorm; bestIdx2 = idx; bestIdx3 = -2; bestIdx4 = -2;
                        }
					}
                    if(bestIdx3==-2) continue;

                    float ag1 = mapmax1[(int)(xf*0.5f+0.25f) + (int)(yf*0.5f+0.25f)*w1];
					if(ag1 > pixelTH1*thFactor)
					{
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir3)));
						if(!setting_selectDirectionDistribution) dirNorm = ag1;

						if(dirNorm > bestVal3)
						{ bestVal3 = dirNorm; bestIdx3 = idx; bestIdx4 = -2;}
					}
					if(bestIdx4==-2) continue;

					float ag2 = mapmax2[(int)(xf*0.25f+0.125) + (int)(yf*0.25f+0.125)*w2];
					if(ag2 > pixelTH2*thFactor)
					{
						Vec2f ag0d = map0[idx].tail<2>();
						float dirNorm = fabsf((float)(ag0d.dot(dir4)));
						if(!setting_selectDirectionDistribution) dirNorm = ag2;

						if(dirNorm > bestVal4)
						{ bestVal4 = dirNorm; bestIdx4 = idx; }
					}                    
                } 

				if(bestIdx2>0)
				{
					map_out[bestIdx2] = 1;
					bestVal3 = 1e10;
					n2++;
				}                             
            }

			if(bestIdx3>0)
			{
				map_out[bestIdx3] = 2;
				bestVal4 = 1e10;
				n3++;
			}            
        }

		if(bestIdx4>0)
		{
			map_out[bestIdx4] = 4;
			n4++;
		}        
    }

    return Eigen::Vector3i(n2,n3,n4);
}