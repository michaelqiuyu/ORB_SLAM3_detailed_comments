/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#include "G2oTypes.h"
#include "ImuTypes.h"
#include "Converter.h"
namespace ORB_SLAM3
{

/** 
 * @brief 储存关键帧位姿相关的信息，用于优化
 * @param pKF 关键帧
 */
ImuCamPose::ImuCamPose(KeyFrame *pKF):its(0)
{
    // Load IMU pose
    twb = pKF->GetImuPosition().cast<double>();
    Rwb = pKF->GetImuRotation().cast<double>();

    // Load camera poses
    int num_cams;
    if(pKF->mpCamera2)
        num_cams=2;
    else
        num_cams=1;

    tcw.resize(num_cams);
    Rcw.resize(num_cams);
    tcb.resize(num_cams);
    Rcb.resize(num_cams);
    Rbc.resize(num_cams);
    tbc.resize(num_cams);
    pCamera.resize(num_cams);

    // Left camera
    tcw[0] = pKF->GetTranslation().cast<double>();
    Rcw[0] = pKF->GetRotation().cast<double>();
    tcb[0] = pKF->mImuCalib.mTcb.translation().cast<double>();
    Rcb[0] = pKF->mImuCalib.mTcb.rotationMatrix().cast<double>();
    Rbc[0] = Rcb[0].transpose();
    tbc[0] = pKF->mImuCalib.mTbc.translation().cast<double>();
    pCamera[0] = pKF->mpCamera;
    bf = pKF->mbf;

    if(num_cams>1)
    {
        // 双目情形，并一直left→right的变换
        Eigen::Matrix4d Trl = pKF->GetRelativePoseTrl().matrix().cast<double>();
        Rcw[1] = Trl.block<3,3>(0,0) * Rcw[0];
        tcw[1] = Trl.block<3,3>(0,0) * tcw[0] + Trl.block<3,1>(0,3);
        tcb[1] = Trl.block<3,3>(0,0) * tcb[0] + Trl.block<3,1>(0,3);
        Rcb[1] = Trl.block<3,3>(0,0) * Rcb[0];
        Rbc[1] = Rcb[1].transpose();
        tbc[1] = -Rbc[1] * tcb[1];
        pCamera[1] = pKF->mpCamera2;
    }

    // For posegraph 4DoF
    Rwb0 = Rwb;
    DR.setIdentity();
}

/** 
 * @brief 储存普通帧位姿相关的信息，用于优化
 * @param pF 普通帧
 */
ImuCamPose::ImuCamPose(Frame *pF):its(0)
{
    // Load IMU pose
    twb = pF->GetImuPosition().cast<double>();
    Rwb = pF->GetImuRotation().cast<double>();

    // Load camera poses
    int num_cams;
    if(pF->mpCamera2)
        num_cams=2;
    else
        num_cams=1;

    tcw.resize(num_cams);
    Rcw.resize(num_cams);
    tcb.resize(num_cams);
    Rcb.resize(num_cams);
    Rbc.resize(num_cams);
    tbc.resize(num_cams);
    pCamera.resize(num_cams);

    // Left camera
    tcw[0] = pF->GetPose().translation().cast<double>();
    Rcw[0] = pF->GetPose().rotationMatrix().cast<double>();
    tcb[0] = pF->mImuCalib.mTcb.translation().cast<double>();
    Rcb[0] = pF->mImuCalib.mTcb.rotationMatrix().cast<double>();
    Rbc[0] = Rcb[0].transpose();
    tbc[0] = pF->mImuCalib.mTbc.translation().cast<double>();
    pCamera[0] = pF->mpCamera;
    bf = pF->mbf;

    if(num_cams>1)
    {
        Eigen::Matrix4d Trl = pF->GetRelativePoseTrl().matrix().cast<double>();
        Rcw[1] = Trl.block<3,3>(0,0) * Rcw[0];
        tcw[1] = Trl.block<3,3>(0,0) * tcw[0] + Trl.block<3,1>(0,3);
        tcb[1] = Trl.block<3,3>(0,0) * tcb[0] + Trl.block<3,1>(0,3);
        Rcb[1] = Trl.block<3,3>(0,0) * Rcb[0];
        Rbc[1] = Rcb[1].transpose();
        tbc[1] = -Rbc[1] * tcb[1];
        pCamera[1] = pF->mpCamera2;
    }

    // For posegraph 4DoF
    Rwb0 = Rwb;
    DR.setIdentity();
}

/** 
 * @brief 储存位姿相关的信息，用于优化
 */
ImuCamPose::ImuCamPose(Eigen::Matrix3d &_Rwc, Eigen::Vector3d &_twc, KeyFrame* pKF): its(0)
{
    // notes: 使用的是传入的R，t，而不是pKF的位姿
    // This is only for posegrpah, we do not care about multicamera
    tcw.resize(1);
    Rcw.resize(1);
    tcb.resize(1);
    Rcb.resize(1);
    Rbc.resize(1);
    tbc.resize(1);
    pCamera.resize(1);

    tcb[0] = pKF->mImuCalib.mTcb.translation().cast<double>();
    Rcb[0] = pKF->mImuCalib.mTcb.rotationMatrix().cast<double>();
    Rbc[0] = Rcb[0].transpose();
    tbc[0] = pKF->mImuCalib.mTbc.translation().cast<double>();
    twb = _Rwc * tcb[0] + _twc;
    Rwb = _Rwc * Rcb[0];
    Rcw[0] = _Rwc.transpose();
    tcw[0] = -Rcw[0] * _twc;
    pCamera[0] = pKF->mpCamera;
    bf = pKF->mbf;

    // For posegraph 4DoF
    Rwb0 = Rwb;
    DR.setIdentity();
}

/** 
 * @brief 设置相关数据
 */
void ImuCamPose::SetParam(
    const std::vector<Eigen::Matrix3d> &_Rcw, const std::vector<Eigen::Vector3d> &_tcw,
    const std::vector<Eigen::Matrix3d> &_Rbc, const std::vector<Eigen::Vector3d> &_tbc, const double &_bf)
{
    Rbc = _Rbc;
    tbc = _tbc;
    Rcw = _Rcw;
    tcw = _tcw;
    const int num_cams = Rbc.size();
    Rcb.resize(num_cams);
    tcb.resize(num_cams);

    for(int i=0; i<tcb.size(); i++)
    {
        Rcb[i] = Rbc[i].transpose();
        tcb[i] = -Rcb[i]*tbc[i];
    }
    // notes: 明确哪些是Eigen，哪些是vector<Eigen>
    Rwb = Rcw[0].transpose()*Rcb[0];
    twb = Rcw[0].transpose()*(tcb[0]-tcw[0]);

    bf = _bf;
}

/** 
 * @brief 单目投影
 */
Eigen::Vector2d ImuCamPose::Project(const Eigen::Vector3d &Xw, int cam_idx) const
{
    Eigen::Vector3d Xc = Rcw[cam_idx] * Xw + tcw[cam_idx];

    return pCamera[cam_idx]->project(Xc);
}

/** 
 * @brief 双目投影，都是0
 * @return u v u`
 */
Eigen::Vector3d ImuCamPose::ProjectStereo(const Eigen::Vector3d &Xw, int cam_idx) const
{
    Eigen::Vector3d Pc = Rcw[cam_idx] * Xw + tcw[cam_idx];
    Eigen::Vector3d pc;
    double invZ = 1/Pc(2);
    pc.head(2) = pCamera[cam_idx]->project(Pc);
    pc(2) = pc(0) - bf*invZ;  // 计算的是ur
    return pc;  // ul vl ur
}

/** 
 * @brief 判断深度值是否有效
 */
bool ImuCamPose::isDepthPositive(const Eigen::Vector3d &Xw, int cam_idx) const
{
    return (Rcw[cam_idx].row(2) * Xw + tcw[cam_idx](2)) > 0.0;  // 相机系下的深度是否为正
}

/** 
 * @brief 优化算出更新值，更新到状态中
 * @param pu 更新值
 */
void ImuCamPose::Update(const double *pu)
{
    Eigen::Vector3d ur, ut;
    ur << pu[0], pu[1], pu[2];
    ut << pu[3], pu[4], pu[5];

    // Update body pose
    // 查看邱晓晨的文档
    twb += Rwb * ut;
    Rwb = Rwb * ExpSO3(ur);

    // Normalize rotation after 5 updates
    its++;
    if(its>=3)
    {
        NormalizeRotation(Rwb);
        its=0;
    }

    // Update camera poses
    const Eigen::Matrix3d Rbw = Rwb.transpose();
    const Eigen::Vector3d tbw = -Rbw * twb;

    for(int i=0; i<pCamera.size(); i++)
    {
        Rcw[i] = Rcb[i] * Rbw;
        tcw[i] = Rcb[i] * tbw + tcb[i];
    }

}

// 更新世界坐标系
void ImuCamPose::UpdateW(const double *pu)
{
    Eigen::Vector3d ur, ut;
    ur << pu[0], pu[1], pu[2];
    ut << pu[3], pu[4], pu[5];


    const Eigen::Matrix3d dR = ExpSO3(ur);
    DR = dR * DR;
    Rwb = DR * Rwb0;
    // Update body pose
    twb += ut;

    // Normalize rotation after 5 updates
    its++;
    if(its>=5)
    {
        // 对于轴角方向为(0, 0, 1)，其对应的旋转矩阵的DR(0,2)、DR(0,2)、DR(0,2)、DR(2,1)均为0
        DR(0,2) = 0.0;
        DR(0,2) = 0.0;
        DR(0,2) = 0.0;
        DR(2,1) = 0.0;
        NormalizeRotation(DR);
        its = 0;
        // 只要修正DR就可以修正Rwb了，因为Rwb0不变
    }

    // Update camera pose
    const Eigen::Matrix3d Rbw = Rwb.transpose();
    const Eigen::Vector3d tbw = -Rbw * twb;

    for(int i=0; i<pCamera.size(); i++)
    {
        Rcw[i] = Rcb[i] * Rbw;
        tcw[i] = Rcb[i] * tbw+tcb[i];
    }
}

// 关于逆深度的，暂未使用
InvDepthPoint::InvDepthPoint(double _rho, double _u, double _v, KeyFrame* pHostKF): u(_u), v(_v), rho(_rho),
    fx(pHostKF->fx), fy(pHostKF->fy), cx(pHostKF->cx), cy(pHostKF->cy), bf(pHostKF->mbf)
{
}

void InvDepthPoint::Update(const double *pu)
{
    rho += *pu;
}

// xc's todo: 在哪里使用了这两个函数
/** 
 * @brief 写入状态量
 */
bool VertexPose::read(std::istream& is)
{
    std::vector<Eigen::Matrix<double,3,3> > Rcw;
    std::vector<Eigen::Matrix<double,3,1> > tcw;
    std::vector<Eigen::Matrix<double,3,3> > Rbc;
    std::vector<Eigen::Matrix<double,3,1> > tbc;

    const int num_cams = _estimate.Rbc.size();
    for(int idx = 0; idx<num_cams; idx++)
    {
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++)
                is >> Rcw[idx](i,j);
        }
        for (int i=0; i<3; i++){
            is >> tcw[idx](i);
        }

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++)
                is >> Rbc[idx](i,j);
        }
        for (int i=0; i<3; i++){
            is >> tbc[idx](i);
        }

        float nextParam;
        for(size_t i = 0; i < _estimate.pCamera[idx]->size(); i++){
            is >> nextParam;
            _estimate.pCamera[idx]->setParameter(nextParam,i);
        }
    }

    double bf;
    is >> bf;
    _estimate.SetParam(Rcw,tcw,Rbc,tbc,bf);
    updateCache();
    
    return true;
}

/** 
 * @brief 读出状态量
 */
bool VertexPose::write(std::ostream& os) const
{
    std::vector<Eigen::Matrix<double,3,3> > Rcw = _estimate.Rcw;
    std::vector<Eigen::Matrix<double,3,1> > tcw = _estimate.tcw;

    std::vector<Eigen::Matrix<double,3,3> > Rbc = _estimate.Rbc;
    std::vector<Eigen::Matrix<double,3,1> > tbc = _estimate.tbc;

    const int num_cams = tcw.size();

    for(int idx = 0; idx<num_cams; idx++)
    {
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++)
                os << Rcw[idx](i,j) << " ";
        }
        for (int i=0; i<3; i++){
            os << tcw[idx](i) << " ";
        }

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++)
                os << Rbc[idx](i,j) << " ";
        }
        for (int i=0; i<3; i++){
            os << tbc[idx](i) << " ";
        }

        for(size_t i = 0; i < _estimate.pCamera[idx]->size(); i++){
            os << _estimate.pCamera[idx]->getParameter(i) << " ";
        }
    }

    os << _estimate.bf << " ";

    return os.good();
}

/** 
 * @brief 单目视觉边计算雅克比
 * _jacobianOplusXi对应着_vertices[0] 也就是误差对于三维点的雅克比
 * _jacobianOplusXj就对应着位姿
 */
void EdgeMono::linearizeOplus()
{
    const VertexPose* VPose = static_cast<const VertexPose*>(_vertices[1]);
    const g2o::VertexSBAPointXYZ* VPoint = static_cast<const g2o::VertexSBAPointXYZ*>(_vertices[0]);

    const Eigen::Matrix3d &Rcw = VPose->estimate().Rcw[cam_idx];
    const Eigen::Vector3d &tcw = VPose->estimate().tcw[cam_idx];
    const Eigen::Vector3d Xc = Rcw*VPoint->estimate() + tcw;
    // 这里的Xb应该也可以直接由Tbw计算得到
    const Eigen::Vector3d Xb = VPose->estimate().Rbc[cam_idx]*Xc+VPose->estimate().tbc[cam_idx];
    const Eigen::Matrix3d &Rcb = VPose->estimate().Rcb[cam_idx];

    /**
     * notes:
     *      1. 优化的是Rbw，因此投影过程为K * Tcb * Tbw * Pw
     *      2. 令p1 = Tcb * Tbw * Pw = Tcw * Pw
     *      3. 令p2 = Tbw * Pw
     *      4. 根据链式求导即可得到下面的结果
     */
    const Eigen::Matrix<double,2,3> proj_jac = VPose->estimate().pCamera[cam_idx]->projectJac(Xc);
    _jacobianOplusXi = -proj_jac * Rcw;

    Eigen::Matrix<double,3,6> SE3deriv;
    double x = Xb(0);
    double y = Xb(1);
    double z = Xb(2);

    SE3deriv <<  0.0,   z,  -y, 1.0, 0.0, 0.0,
                 -z , 0.0,   x, 0.0, 1.0, 0.0,
                  y ,  -x, 0.0, 0.0, 0.0, 1.0;

    _jacobianOplusXj = proj_jac * Rcb * SE3deriv; // TODO optimize this product
}

/** 
 * @brief 单目视觉纯位姿边计算雅克比
 * _jacobianOplusXi对应着_vertices[0] 也就是误差对于位姿的雅克比
 */
void EdgeMonoOnlyPose::linearizeOplus()
{
    const VertexPose* VPose = static_cast<const VertexPose*>(_vertices[0]);

    const Eigen::Matrix3d &Rcw = VPose->estimate().Rcw[cam_idx];
    const Eigen::Vector3d &tcw = VPose->estimate().tcw[cam_idx];
    const Eigen::Vector3d Xc = Rcw*Xw + tcw;
    const Eigen::Vector3d Xb = VPose->estimate().Rbc[cam_idx]*Xc+VPose->estimate().tbc[cam_idx];
    const Eigen::Matrix3d &Rcb = VPose->estimate().Rcb[cam_idx];

    Eigen::Matrix<double,2,3> proj_jac = VPose->estimate().pCamera[cam_idx]->projectJac(Xc);

    Eigen::Matrix<double,3,6> SE3deriv;
    double x = Xb(0);
    double y = Xb(1);
    double z = Xb(2);
    SE3deriv << 0.0,   z,   -y, 1.0, 0.0, 0.0,
                -z , 0.0,    x, 0.0, 1.0, 0.0,
                 y , - x,  0.0, 0.0, 0.0, 1.0;
    _jacobianOplusXi = proj_jac * Rcb * SE3deriv; // symbol different becasue of update mode
}

/** 
 * @brief 双目视觉边计算雅克比，多了一维误差
 * _jacobianOplusXi对应着_vertices[0] 也就是误差对于三维点的雅克比
 * _jacobianOplusXj就对应着位姿
 */
void EdgeStereo::linearizeOplus()
{
    const VertexPose* VPose = static_cast<const VertexPose*>(_vertices[1]);
    const g2o::VertexSBAPointXYZ* VPoint = static_cast<const g2o::VertexSBAPointXYZ*>(_vertices[0]);

    const Eigen::Matrix3d &Rcw = VPose->estimate().Rcw[cam_idx];
    const Eigen::Vector3d &tcw = VPose->estimate().tcw[cam_idx];
    const Eigen::Vector3d Xc = Rcw*VPoint->estimate() + tcw;
    const Eigen::Vector3d Xb = VPose->estimate().Rbc[cam_idx]*Xc+VPose->estimate().tbc[cam_idx];
    const Eigen::Matrix3d &Rcb = VPose->estimate().Rcb[cam_idx];
    const double bf = VPose->estimate().bf;
    const double inv_z2 = 1.0/(Xc(2)*Xc(2));

    /**
     * ur = ul - fb / z
     * 根据上面的公式求解ur对x, y, z的导数
     */
    Eigen::Matrix<double,3,3> proj_jac;
    proj_jac.block<2,3>(0,0) = VPose->estimate().pCamera[cam_idx]->projectJac(Xc);
    proj_jac.block<1,3>(2,0) = proj_jac.block<1,3>(0,0);
    proj_jac(2,2) += bf*inv_z2;

    _jacobianOplusXi = -proj_jac * Rcw;

    Eigen::Matrix<double,3,6> SE3deriv;
    double x = Xb(0);
    double y = Xb(1);
    double z = Xb(2);

    SE3deriv << 0.0, z,   -y, 1.0, 0.0, 0.0,
            -z , 0.0, x, 0.0, 1.0, 0.0,
            y ,  -x , 0.0, 0.0, 0.0, 1.0;

    _jacobianOplusXj = proj_jac * Rcb * SE3deriv;
}

/** 
 * @brief 双目视觉纯位姿边计算雅克比，多了一维误差
 * _jacobianOplusXi对应着_vertices[0] 也就是误差对于位姿的雅克比
 */
void EdgeStereoOnlyPose::linearizeOplus()
{
    const VertexPose* VPose = static_cast<const VertexPose*>(_vertices[0]);

    const Eigen::Matrix3d &Rcw = VPose->estimate().Rcw[cam_idx];
    const Eigen::Vector3d &tcw = VPose->estimate().tcw[cam_idx];
    const Eigen::Vector3d Xc = Rcw*Xw + tcw;
    const Eigen::Vector3d Xb = VPose->estimate().Rbc[cam_idx]*Xc+VPose->estimate().tbc[cam_idx];
    const Eigen::Matrix3d &Rcb = VPose->estimate().Rcb[cam_idx];
    const double bf = VPose->estimate().bf;
    const double inv_z2 = 1.0/(Xc(2)*Xc(2));

    Eigen::Matrix<double,3,3> proj_jac;
    proj_jac.block<2,3>(0,0) = VPose->estimate().pCamera[cam_idx]->projectJac(Xc);
    proj_jac.block<1,3>(2,0) = proj_jac.block<1,3>(0,0);
    proj_jac(2,2) += bf*inv_z2;

    Eigen::Matrix<double,3,6> SE3deriv;
    double x = Xb(0);
    double y = Xb(1);
    double z = Xb(2);
    SE3deriv << 0.0, z,   -y, 1.0, 0.0, 0.0,
            -z , 0.0, x, 0.0, 1.0, 0.0,
            y ,  -x , 0.0, 0.0, 0.0, 1.0;
    _jacobianOplusXi = proj_jac * Rcb * SE3deriv;
}

VertexVelocity::VertexVelocity(KeyFrame* pKF)
{
    setEstimate(pKF->GetVelocity().cast<double>());
}

VertexVelocity::VertexVelocity(Frame* pF)
{
    setEstimate(pF->GetVelocity().cast<double>());
}

VertexGyroBias::VertexGyroBias(KeyFrame *pKF)
{
    setEstimate(pKF->GetGyroBias().cast<double>());
}

VertexGyroBias::VertexGyroBias(Frame *pF)
{
    Eigen::Vector3d bg;
    bg << pF->mImuBias.bwx, pF->mImuBias.bwy,pF->mImuBias.bwz;
    setEstimate(bg);
}

VertexAccBias::VertexAccBias(KeyFrame *pKF)
{
    setEstimate(pKF->GetAccBias().cast<double>());
}

VertexAccBias::VertexAccBias(Frame *pF)
{
    Eigen::Vector3d ba;
    ba << pF->mImuBias.bax, pF->mImuBias.bay,pF->mImuBias.baz;
    setEstimate(ba);
}


/** 
 * @brief 局部地图中imu的局部地图优化（此时已经初始化完毕不需要再优化重力方向与尺度）
 * @param pInt 预积分相关内容
 */
EdgeInertial::EdgeInertial(IMU::Preintegrated *pInt):JRg(pInt->JRg.cast<double>()),
    JVg(pInt->JVg.cast<double>()), JPg(pInt->JPg.cast<double>()), JVa(pInt->JVa.cast<double>()),
    JPa(pInt->JPa.cast<double>()), mpInt(pInt), dt(pInt->dT)
{
    // 准备工作，把预积分类里面的值先取出来，包含信息的是两帧之间n多个imu信息预积分来的
    // This edge links 6 vertices
    // 6元边
    resize(6);
    // 1. 定义重力
    g << 0, 0, -IMU::GRAVITY_VALUE;

    // 2. 读取协方差矩阵的前9*9部分的逆矩阵，该部分表示的是预积分测量噪声的协方差矩阵
    /**
     * notes:
     *      1. 这里C.block<9, 9>不一定是非奇异矩阵，因此这里的操作应该修改为对C.block<9, 9>做特征值分解，并利用特征值分解求逆
     *      2. 也就是对A求特征值分解，A = p * diag * p.t→A.inv = p * diag.inv * p，且diag(i) = 0时，令diag.inv(i) = 0，这样做的理由如下
     *
     * 为了说明这样操作的理由，不妨举一个特殊的例子，info = diag(a1, a2, ..., an)，ai表示方差，如果ai≠0，那么
     * J = (e1.t, e2.t, ..., en.t).t，J.t * info * J = ∑ei.t * ei / ai
     * ai等于0，也就是方差为0，也就是没有波动，此时不考虑误差项ei.t * ei，因此直接将1/ai令为0即可
     */
    Matrix9d Info = pInt->C.block<9,9>(0,0).cast<double>().inverse();
    // 3. 强制让其成为对角矩阵
    // 理论上Info就是一个对称矩阵，由于浮点数误差，可能变成非对称矩阵了
    Info = (Info+Info.transpose())/2;
    // 4. 让特征值很小的时候置为0，再重新计算信息矩阵（暂不知这么操作的目的是什么，先搞清楚操作流程吧）
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,9,9> > es(Info);
    Eigen::Matrix<double,9,1> eigs = es.eigenvalues();
    for(int i=0;i<9;i++)
        if(eigs[i]<1e-12)
            // 太小，说明矩阵退化了，这里的很小的浮点数可能是浮点数运算误差带来的，直接变为0即可，也可以避免某些无用的运算
            eigs[i]=0;
    // asDiagonal 生成对角矩阵
    Info = es.eigenvectors()*eigs.asDiagonal()*es.eigenvectors().transpose();
    setInformation(Info);
}



/** 
 * @brief 计算误差
 */
void EdgeInertial::computeError()
{
    // TODO Maybe Reintegrate inertial measurments when difference between linearization point and current estimate is too big
    const VertexPose* VP1 = static_cast<const VertexPose*>(_vertices[0]);           //位姿Ti
    const VertexVelocity* VV1= static_cast<const VertexVelocity*>(_vertices[1]);    //速度vi
    const VertexGyroBias* VG1= static_cast<const VertexGyroBias*>(_vertices[2]);    //零偏Bgi
    const VertexAccBias* VA1= static_cast<const VertexAccBias*>(_vertices[3]);      //零偏Bai
    const VertexPose* VP2 = static_cast<const VertexPose*>(_vertices[4]);           //位姿Tj
    const VertexVelocity* VV2 = static_cast<const VertexVelocity*>(_vertices[5]);   //速度vj
    const IMU::Bias b1(VA1->estimate()[0],VA1->estimate()[1],VA1->estimate()[2],VG1->estimate()[0],VG1->estimate()[1],VG1->estimate()[2]);
    const Eigen::Matrix3d dR = mpInt->GetDeltaRotation(b1).cast<double>();
    const Eigen::Vector3d dV = mpInt->GetDeltaVelocity(b1).cast<double>();
    const Eigen::Vector3d dP = mpInt->GetDeltaPosition(b1).cast<double>();

    const Eigen::Vector3d er = LogSO3(dR.transpose()*VP1->estimate().Rwb.transpose()*VP2->estimate().Rwb);
    const Eigen::Vector3d ev = VP1->estimate().Rwb.transpose()*(VV2->estimate() - VV1->estimate() - g*dt) - dV;
    const Eigen::Vector3d ep = VP1->estimate().Rwb.transpose()*(VP2->estimate().twb - VP1->estimate().twb
                                                               - VV1->estimate()*dt - g*dt*dt/2) - dP;

    _error << er, ev, ep;
}

// 计算雅克比矩阵
void EdgeInertial::linearizeOplus()
{
    const VertexPose* VP1 = static_cast<const VertexPose*>(_vertices[0]);
    const VertexVelocity* VV1= static_cast<const VertexVelocity*>(_vertices[1]);
    const VertexGyroBias* VG1= static_cast<const VertexGyroBias*>(_vertices[2]);
    const VertexAccBias* VA1= static_cast<const VertexAccBias*>(_vertices[3]);
    const VertexPose* VP2 = static_cast<const VertexPose*>(_vertices[4]);
    const VertexVelocity* VV2= static_cast<const VertexVelocity*>(_vertices[5]);
    const IMU::Bias b1(VA1->estimate()[0],VA1->estimate()[1],VA1->estimate()[2],VG1->estimate()[0],VG1->estimate()[1],VG1->estimate()[2]);
    const IMU::Bias db = mpInt->GetDeltaBias(b1);
    Eigen::Vector3d dbg;
    dbg << db.bwx, db.bwy, db.bwz;

    const Eigen::Matrix3d Rwb1 = VP1->estimate().Rwb;  // Ri
    const Eigen::Matrix3d Rbw1 = Rwb1.transpose();     // Ri.t()
    const Eigen::Matrix3d Rwb2 = VP2->estimate().Rwb;  // Rj

    const Eigen::Matrix3d dR = mpInt->GetDeltaRotation(b1).cast<double>();
    const Eigen::Matrix3d eR = dR.transpose()*Rbw1*Rwb2;        // r△Rij
    const Eigen::Vector3d er = LogSO3(eR);                      // r△φij
    const Eigen::Matrix3d invJr = InverseRightJacobianSO3(er);  // Jr^-1(log(△Rij))

    // 这里可以直接参考邱笑晨的文档，公式完全一样
    // 就很神奇，_jacobianOplus个数等于边的个数，里面的大小等于观测值维度（也就是残差）× 每个节点待优化值的维度
    // Jacobians wrt Pose 1
    // _jacobianOplus[0] 9*6矩阵 总体来说就是三个残差分别对pose1的旋转与平移（p）求导
    _jacobianOplus[0].setZero();
    // rotation
    // (0,0)起点的3*3块表示旋转残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(0,0) = -invJr*Rwb2.transpose()*Rwb1; // OK
    // (3,0)起点的3*3块表示速度残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(3,0) = Sophus::SO3d::hat(Rbw1*(VV2->estimate() - VV1->estimate() - g*dt)); // OK
    // (6,0)起点的3*3块表示位置残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(6,0) = Sophus::SO3d::hat(Rbw1*(VP2->estimate().twb - VP1->estimate().twb
                                                   - VV1->estimate()*dt - 0.5*g*dt*dt)); // OK
    // translation
    // (6,3)起点的3*3块表示位置残差对pose1的位置求导
    _jacobianOplus[0].block<3,3>(6,3) = -Eigen::Matrix3d::Identity(); // OK

    // Jacobians wrt Velocity 1
    // _jacobianOplus[1] 9*3矩阵 总体来说就是三个残差分别对pose1的速度求导
    _jacobianOplus[1].setZero();
    _jacobianOplus[1].block<3,3>(3,0) = -Rbw1; // OK
    _jacobianOplus[1].block<3,3>(6,0) = -Rbw1*dt; // OK

    // Jacobians wrt Gyro 1
    // _jacobianOplus[2] 9*3矩阵 总体来说就是三个残差分别对陀螺仪偏置的速度求导
    _jacobianOplus[2].setZero();
    _jacobianOplus[2].block<3,3>(0,0) = -invJr*eR.transpose()*RightJacobianSO3(JRg*dbg)*JRg; // OK
    _jacobianOplus[2].block<3,3>(3,0) = -JVg; // OK
    _jacobianOplus[2].block<3,3>(6,0) = -JPg; // OK

    // Jacobians wrt Accelerometer 1
    // _jacobianOplus[3] 9*3矩阵 总体来说就是三个残差分别对加速度计偏置的速度求导
    _jacobianOplus[3].setZero();
    _jacobianOplus[3].block<3,3>(3,0) = -JVa; // OK
    _jacobianOplus[3].block<3,3>(6,0) = -JPa; // OK

    // Jacobians wrt Pose 2
    // _jacobianOplus[4] 9*6矩阵 总体来说就是三个残差分别对pose2的旋转与平移（p）求导
    _jacobianOplus[4].setZero();
    // rotation
    _jacobianOplus[4].block<3,3>(0,0) = invJr; // OK
    // translation
    _jacobianOplus[4].block<3,3>(6,3) = Rbw1*Rwb2; // OK

    // Jacobians wrt Velocity 2
    // _jacobianOplus[5] 9*3矩阵 总体来说就是三个残差分别对pose2的速度求导
    _jacobianOplus[5].setZero();
    _jacobianOplus[5].block<3,3>(3,0) = Rbw1; // OK
}

// localmapping中imu初始化所用的边，除了正常的几个优化变量外还优化了重力方向与尺度
EdgeInertialGS::EdgeInertialGS(IMU::Preintegrated *pInt):JRg(pInt->JRg.cast<double>()),
    JVg(pInt->JVg.cast<double>()), JPg(pInt->JPg.cast<double>()), JVa(pInt->JVa.cast<double>()),
    JPa(pInt->JPa.cast<double>()), mpInt(pInt), dt(pInt->dT)
{
    // 准备工作，把预积分类里面的值先取出来，包含信息的是两帧之间n多个imu信息预积分来的
    // This edge links 8 vertices
    // 8元边
    resize(8);  // 多了尺度和重力
    // 1. 定义重力
    gI << 0, 0, -IMU::GRAVITY_VALUE;

    // 2. 读取协方差矩阵的前9*9部分的逆矩阵，该部分表示的是预积分测量噪声的协方差矩阵
    Matrix9d Info = pInt->C.block<9,9>(0,0).cast<double>().inverse();
    // 3. 强制让其成为对角矩阵
    Info = (Info+Info.transpose())/2;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,9,9> > es(Info);
    // 4. 让特征值很小的时候置为0，再重新计算信息矩阵（暂不知这么操作的目的是什么，先搞清楚操作流程吧）
    Eigen::Matrix<double,9,1> eigs = es.eigenvalues();
    for(int i=0;i<9;i++)
        if(eigs[i]<1e-12)
            eigs[i]=0;
    // asDiagonal 生成对角矩阵
    Info = es.eigenvectors()*eigs.asDiagonal()*es.eigenvectors().transpose();
    setInformation(Info);
}


// 计算误差
void EdgeInertialGS::computeError()
{
    // TODO Maybe Reintegrate inertial measurments when difference between linearization point and current estimate is too big
    const VertexPose* VP1 = static_cast<const VertexPose*>(_vertices[0]);
    const VertexVelocity* VV1= static_cast<const VertexVelocity*>(_vertices[1]);
    const VertexGyroBias* VG= static_cast<const VertexGyroBias*>(_vertices[2]);
    const VertexAccBias* VA= static_cast<const VertexAccBias*>(_vertices[3]);
    const VertexPose* VP2 = static_cast<const VertexPose*>(_vertices[4]);
    const VertexVelocity* VV2 = static_cast<const VertexVelocity*>(_vertices[5]);
    const VertexGDir* VGDir = static_cast<const VertexGDir*>(_vertices[6]);
    const VertexScale* VS = static_cast<const VertexScale*>(_vertices[7]);
    const IMU::Bias b(VA->estimate()[0],VA->estimate()[1],VA->estimate()[2],VG->estimate()[0],VG->estimate()[1],VG->estimate()[2]);
    g = VGDir->estimate().Rwg*gI;
    const double s = VS->estimate();
    const Eigen::Matrix3d dR = mpInt->GetDeltaRotation(b).cast<double>();
    const Eigen::Vector3d dV = mpInt->GetDeltaVelocity(b).cast<double>();
    const Eigen::Vector3d dP = mpInt->GetDeltaPosition(b).cast<double>();

    // 计算残差。广义上讲都是真实值 = 残差 + imu，旋转为imu*残差=真实值
    // dR.transpose() 为imu预积分的值，VP1->estimate().Rwb.transpose() * VP2->estimate().Rwb 为相机的Rwc在乘上相机与imu的标定外参矩阵
    const Eigen::Vector3d er = LogSO3(dR.transpose()*VP1->estimate().Rwb.transpose()*VP2->estimate().Rwb);
    const Eigen::Vector3d ev = VP1->estimate().Rwb.transpose()*(s*(VV2->estimate() - VV1->estimate()) - g*dt) - dV;
    const Eigen::Vector3d ep = VP1->estimate().Rwb.transpose()*(s*(VP2->estimate().twb - VP1->estimate().twb - VV1->estimate()*dt) - g*dt*dt/2) - dP;

    _error << er, ev, ep;
}

// 计算雅克比矩阵
// 对重力和尺度方面的导数参考Inertial-Only ptimizatin for Visual-Inertial Initialization的公式6-11
void EdgeInertialGS::linearizeOplus()
{
    const VertexPose* VP1 = static_cast<const VertexPose*>(_vertices[0]);
    const VertexVelocity* VV1= static_cast<const VertexVelocity*>(_vertices[1]);
    const VertexGyroBias* VG= static_cast<const VertexGyroBias*>(_vertices[2]);
    const VertexAccBias* VA= static_cast<const VertexAccBias*>(_vertices[3]);
    const VertexPose* VP2 = static_cast<const VertexPose*>(_vertices[4]);
    const VertexVelocity* VV2 = static_cast<const VertexVelocity*>(_vertices[5]);
    const VertexGDir* VGDir = static_cast<const VertexGDir*>(_vertices[6]);
    const VertexScale* VS = static_cast<const VertexScale*>(_vertices[7]);
    // 1. 获取偏置的该变量，因为要对这个东西求导
    const IMU::Bias b(VA->estimate()[0],VA->estimate()[1],VA->estimate()[2],VG->estimate()[0],VG->estimate()[1],VG->estimate()[2]);
    const IMU::Bias db = mpInt->GetDeltaBias(b);

    // 陀螺仪的偏置改变量
    Eigen::Vector3d dbg;
    dbg << db.bwx, db.bwy, db.bwz;

    const Eigen::Matrix3d Rwb1 = VP1->estimate().Rwb;   // Ri
    const Eigen::Matrix3d Rbw1 = Rwb1.transpose();      // Ri.t()
    const Eigen::Matrix3d Rwb2 = VP2->estimate().Rwb;   // Rj
    const Eigen::Matrix3d Rwg = VGDir->estimate().Rwg;  // Rwg
    Eigen::MatrixXd Gm = Eigen::MatrixXd::Zero(3,2);
    Gm(0,1) = -IMU::GRAVITY_VALUE;
    Gm(1,0) = IMU::GRAVITY_VALUE;
    const double s = VS->estimate();
    const Eigen::MatrixXd dGdTheta = Rwg*Gm;
    const Eigen::Matrix3d dR = mpInt->GetDeltaRotation(b).cast<double>();
    const Eigen::Matrix3d eR = dR.transpose()*Rbw1*Rwb2;        // r△Rij
    const Eigen::Vector3d er = LogSO3(eR);                      // r△φij
    const Eigen::Matrix3d invJr = InverseRightJacobianSO3(er);  // Jr^-1(log(△Rij))

    // 就很神奇，_jacobianOplus个数等于边的个数，里面的大小等于观测值维度（也就是残差）× 每个节点待优化值的维度
    // Jacobians wrt Pose 1
    // _jacobianOplus[0] 9*6矩阵 总体来说就是三个残差分别对pose1的旋转与平移（p）求导
    _jacobianOplus[0].setZero();
    // rotation
    // (0,0)起点的3*3块表示旋转残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(0,0) = -invJr*Rwb2.transpose()*Rwb1;
    // (3,0)起点的3*3块表示速度残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(3,0) = Sophus::SO3d::hat(Rbw1*(s*(VV2->estimate() - VV1->estimate()) - g*dt));
    // (6,0)起点的3*3块表示位置残差对pose1的旋转求导
    _jacobianOplus[0].block<3,3>(6,0) = Sophus::SO3d::hat(Rbw1*(s*(VP2->estimate().twb - VP1->estimate().twb
                                                   - VV1->estimate()*dt) - 0.5*g*dt*dt));
    // translation
    // (6,3)起点的3*3块表示位置残差对pose1的位置求导
    _jacobianOplus[0].block<3,3>(6,3) = Eigen::DiagonalMatrix<double,3>(-s,-s,-s);

    // Jacobians wrt Velocity 1
    // _jacobianOplus[1] 9*3矩阵 总体来说就是三个残差分别对pose1的速度求导
    _jacobianOplus[1].setZero();
    _jacobianOplus[1].block<3,3>(3,0) = -s*Rbw1;
    _jacobianOplus[1].block<3,3>(6,0) = -s*Rbw1*dt;

    // Jacobians wrt Gyro bias
    // _jacobianOplus[2] 9*3矩阵 总体来说就是三个残差分别对陀螺仪偏置的速度求导
    _jacobianOplus[2].setZero();
    _jacobianOplus[2].block<3,3>(0,0) = -invJr*eR.transpose()*RightJacobianSO3(JRg*dbg)*JRg;
    _jacobianOplus[2].block<3,3>(3,0) = -JVg;
    _jacobianOplus[2].block<3,3>(6,0) = -JPg;

    // Jacobians wrt Accelerometer bias
    // _jacobianOplus[3] 9*3矩阵 总体来说就是三个残差分别对加速度计偏置的速度求导
    _jacobianOplus[3].setZero();
    _jacobianOplus[3].block<3,3>(3,0) = -JVa;
    _jacobianOplus[3].block<3,3>(6,0) = -JPa;

    // Jacobians wrt Pose 2
    // _jacobianOplus[3] 9*6矩阵 总体来说就是三个残差分别对pose2的旋转与平移（p）求导
    _jacobianOplus[4].setZero();
    // rotation
    _jacobianOplus[4].block<3,3>(0,0) = invJr;
    // translation
    _jacobianOplus[4].block<3,3>(6,3) = s*Rbw1*Rwb2;

    // Jacobians wrt Velocity 2
    // _jacobianOplus[3] 9*3矩阵 总体来说就是三个残差分别对pose2的速度求导
    _jacobianOplus[5].setZero();
    _jacobianOplus[5].block<3,3>(3,0) = s*Rbw1;

    // Jacobians wrt Gravity direction
    // _jacobianOplus[3] 9*2矩阵 总体来说就是三个残差分别对重力方向求导
    _jacobianOplus[6].setZero();
    /**
     * 1. 首先先计算3*3的雅克比矩阵，然后去掉最后一列即可
     * 2. 对于B = A*P而言，去掉B的最后一列也就是去掉P的最后一列，这也是GM矩阵的由来
     */
    _jacobianOplus[6].block<3,2>(3,0) = -Rbw1*dGdTheta*dt;
    _jacobianOplus[6].block<3,2>(6,0) = -0.5*Rbw1*dGdTheta*dt*dt;

    // Jacobians wrt scale factor
    // _jacobianOplus[3] 9*1矩阵 总体来说就是三个残差分别对尺度求导
    /**
     * notes:
     *      1. 疑似BUG：根据尺度更新的规则，这里的雅克比计算的应该有问题，应该再乘以一个尺度才能跟尺度更新的规则一致
     */
    _jacobianOplus[7].setZero();
    _jacobianOplus[7].block<3,1>(3,0) = Rbw1*(VV2->estimate()-VV1->estimate());
    _jacobianOplus[7].block<3,1>(6,0) = Rbw1*(VP2->estimate().twb-VP1->estimate().twb-VV1->estimate()*dt);
}

/** 
 * @brief 滑窗边缘化时用的先验边
 */
EdgePriorPoseImu::EdgePriorPoseImu(ConstraintPoseImu *c)
{
    resize(4);
    Rwb = c->Rwb;
    twb = c->twb;
    vwb = c->vwb;
    bg = c->bg;
    ba = c->ba;
    setInformation(c->H);
}

/** 
 * @brief 先验边计算误差
 */
void EdgePriorPoseImu::computeError()
{
    const VertexPose* VP = static_cast<const VertexPose*>(_vertices[0]);
    const VertexVelocity* VV = static_cast<const VertexVelocity*>(_vertices[1]);
    const VertexGyroBias* VG = static_cast<const VertexGyroBias*>(_vertices[2]);
    const VertexAccBias* VA = static_cast<const VertexAccBias*>(_vertices[3]);

    const Eigen::Vector3d er = LogSO3(Rwb.transpose()*VP->estimate().Rwb);
    /**
     * notes:
     *      1. 在论文Visual-Inertial Monocular SLAM with Map Reuse中实际上没有et这个约束
     *      2. 这里实际上表达的是基于已知的Rwb，构建两个坐标系T1、T2，他们的旋转均为Rwb，平移分别是VP->estimate().twb、twb
     *      3. 分别计算世界系的原点在这两个坐标系下的坐标，然后认为二者应该足够接近
     *      4. 这个约束实在有些不伦不类，如果固定旋转为一个已知的Rwb，这个约束实际上就是e = VP->estimate().twb-twb，也就是IMU中心在世界系下位置的变动
     *      5. 我们知道当表示为Rbw的时候，IMU中心在世界系的坐标为-Rbw.t * tbw，可以基于IMU中心在世界系的位置来构建残差
     *      6. 但是这里旋转的表示为Rwb，也就没必要这么做了，不然残差的含义就是某个点在两个坐标系下的坐标应该相同，这不太符合常理，一般应该是两个点在同一个坐标系的坐标应该相同
     *      7. 仅讨论此处的话，因为使用了一个固定的Rwb，因此运行起来不会有什么问题，但最好还是去掉Rwb，直接使用e = VP->estimate().twb-twb
     */
    const Eigen::Vector3d et = Rwb.transpose()*(VP->estimate().twb-twb);
    const Eigen::Vector3d ev = VV->estimate() - vwb;
    const Eigen::Vector3d ebg = VG->estimate() - bg;
    const Eigen::Vector3d eba = VA->estimate() - ba;

    _error << er, et, ev, ebg, eba;
}

/** 
 * @brief 先验边计算雅克比
 */
void EdgePriorPoseImu::linearizeOplus()
{
    const VertexPose* VP = static_cast<const VertexPose*>(_vertices[0]);
    const Eigen::Vector3d er = LogSO3(Rwb.transpose()*VP->estimate().Rwb);
    // 就很神奇，_jacobianOplus个数等于边的个数，里面的大小等于观测值维度（也就是3旋转3平移3速度6偏置）× 每个节点待优化值的维度
    // 源码可读性太差了。。。里面会自动分配矩阵大小，计算改变量时按照对应位置来
    _jacobianOplus[0].setZero();
    // LOG(Rbw*R*EXP(φ)) = LOG(EXP(LOG(Rbw*R) + Jr(-1)*φ)) = LOG(Rbw*R) + Jr(-1)*φ
    _jacobianOplus[0].block<3,3>(0,0) = InverseRightJacobianSO3(er);   // Jr(-1)
    // Rbw*(t + R*δt - twb) = Rbw*(t - twb) + Rbw*R*δt，注意平移的更新方式
    _jacobianOplus[0].block<3,3>(3,3) = Rwb.transpose()*VP->estimate().Rwb;  // Rbw*R
    _jacobianOplus[1].setZero();
    _jacobianOplus[1].block<3,3>(6,0) = Eigen::Matrix3d::Identity();
    _jacobianOplus[2].setZero();
    _jacobianOplus[2].block<3,3>(9,0) = Eigen::Matrix3d::Identity();
    _jacobianOplus[3].setZero();
    _jacobianOplus[3].block<3,3>(12,0) = Eigen::Matrix3d::Identity();
}

void EdgePriorAcc::linearizeOplus()
{
    // Jacobian wrt bias
    // notes: 疑似BUG：雅克比为-I
    _jacobianOplusXi.block<3,3>(0,0) = Eigen::Matrix3d::Identity();

}

void EdgePriorGyro::linearizeOplus()
{
    // Jacobian wrt bias
    // notes: 疑似BUG：雅克比为-I
    _jacobianOplusXi.block<3,3>(0,0) = Eigen::Matrix3d::Identity();

}

// SO3 FUNCTIONS
Eigen::Matrix3d ExpSO3(const Eigen::Vector3d &w)
{
    return ExpSO3(w[0],w[1],w[2]);
}

Eigen::Matrix3d ExpSO3(const double x, const double y, const double z)
{
    const double d2 = x*x+y*y+z*z;
    const double d = sqrt(d2);
    Eigen::Matrix3d W;
    W << 0.0, -z, y,z, 0.0, -x,-y,  x, 0.0;
    if(d<1e-5)
    {
        Eigen::Matrix3d res = Eigen::Matrix3d::Identity() + W +0.5*W*W;
        return NormalizeRotation(res);
    }
    else
    {
        Eigen::Matrix3d res =Eigen::Matrix3d::Identity() + W*sin(d)/d + W*W*(1.0-cos(d))/d2;
        return NormalizeRotation(res);
    }
}

Eigen::Vector3d LogSO3(const Eigen::Matrix3d &R)
{
    const double tr = R(0,0)+R(1,1)+R(2,2);
    Eigen::Vector3d w;
    w << (R(2,1)-R(1,2))/2, (R(0,2)-R(2,0))/2, (R(1,0)-R(0,1))/2;
    const double costheta = (tr-1.0)*0.5f;
    // xc's todo: 这里costheta < -1的时候直接取w的理由是什么？
    if(costheta>1 || costheta<-1)
        return w;
    const double theta = acos(costheta);
    const double s = sin(theta);
    if(fabs(s)<1e-5)
        return w;
    else
        return theta*w/s;
}

Eigen::Matrix3d InverseRightJacobianSO3(const Eigen::Vector3d &v)
{
    return InverseRightJacobianSO3(v[0],v[1],v[2]);
}

Eigen::Matrix3d InverseRightJacobianSO3(const double x, const double y, const double z)
{
    const double d2 = x*x+y*y+z*z;
    const double d = sqrt(d2);

    Eigen::Matrix3d W;
    W << 0.0, -z, y,z, 0.0, -x,-y,  x, 0.0;
    if(d<1e-5)
        return Eigen::Matrix3d::Identity();
    else
        return Eigen::Matrix3d::Identity() + W/2 + W*W*(1.0/d2 - (1.0+cos(d))/(2.0*d*sin(d)));
}

Eigen::Matrix3d RightJacobianSO3(const Eigen::Vector3d &v)
{
    return RightJacobianSO3(v[0],v[1],v[2]);
}

Eigen::Matrix3d RightJacobianSO3(const double x, const double y, const double z)
{
    const double d2 = x*x+y*y+z*z;
    const double d = sqrt(d2);

    Eigen::Matrix3d W;
    W << 0.0, -z, y,z, 0.0, -x,-y,  x, 0.0;
    if(d<1e-5)
    {
        return Eigen::Matrix3d::Identity();
    }
    else
    {
        return Eigen::Matrix3d::Identity() - W*(1.0-cos(d))/d2 + W*W*(d-sin(d))/(d2*d);
    }
}

Eigen::Matrix3d Skew(const Eigen::Vector3d &w)
{
    Eigen::Matrix3d W;
    W << 0.0, -w[2], w[1],w[2], 0.0, -w[0],-w[1],  w[0], 0.0;
    return W;
}

}
