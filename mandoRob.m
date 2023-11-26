%########### DH MDH 机械臂运动学  
%########### 机械臂动力学 M(q)*ddq + C(q,dq)dq + G(q) = u

classdef mandoRob <handle

    properties 
        DH
        MDH
        n
        R
        T
        Tn
        P
        q_sym
        TT
        M
        C
        G
    
    end
    
    methods
        % 构造函数 全为旋转关节
        % 默认使用MDH参数计算运动学
        function  self=mandoRob(dh,varargin)
    
            alpha = dh(1,:);
            a = dh(2,:);
            d = dh(3,:);
            theta = dh(4,:);
            n = size(dh,2);     
            
            self.n=n;
            q_sym = sym('q%d',[n,1],'real');
            self.q_sym=q_sym;
    
            Tn=eye(4);
      
    %% 
            if strcmp(varargin{1},'DH')
                self.DH=dh;
                for i=1:n
                Rot_theta_q=[cos(q_sym(i)) -sin(q_sym(i)) 0 0;sin(q_sym(i)) cos(q_sym(i)) 0 0;0 0 1 0;0 0 0 1];
    %############### DH参数 单个转移矩阵表示为
    %############### A=Rot(theta)*Tran(d)*Tran(a)*Rot(alpha)    
                Rot_theta=[cos(theta(i)) -sin(theta(i)) 0 0;sin(theta(i)) cos(theta(i)) 0 0;0 0 1 0;0 0 0 1];
                Tran_d=[1 0 0 0; 0 1 0 0; 0 0 1 d(i); 0 0 0 1];
                Tran_a=[1 0 0 a(i); 0 1 0 0; 0 0 1 0; 0 0 0 1];
                Rot_alpha=[1 0 0 0; 0 cos(alpha(i)) -sin(alpha(i)) 0; 0 sin(alpha(i)) cos(alpha(i)) 0; 0 0 0 1];
                T(:,:,i)=simplify(Rot_theta_q*Rot_theta*Tran_d*Tran_a*Rot_alpha);         
                R(:,:,i)=T(1:3,1:3,i);
                P(:,i)=T(1:3,4,i);
                Tn=Tn*T(:,:,i);
                TT{i}=Tn;
                end
                self.T=T;
                self.Tn=Tn;
                self.P=P;
                self.R=R;
                self.TT=TT;
            end
            %% 
            if strcmp(varargin{1},'MDH')
                self.MDH=dh;
                for i=1:n
                Rot_theta_q=[cos(q_sym(i)) -sin(q_sym(i)) 0 0;sin(q_sym(i)) cos(q_sym(i)) 0 0;0 0 1 0;0 0 0 1];
    %############### MDH参数 单个转移矩阵表示为
    %############### A=Rot(alpha)*Tran(a)*Tran(d)*Rot(theta)    
                Rot_theta=[cos(theta(i)) -sin(theta(i)) 0 0;sin(theta(i)) cos(theta(i)) 0 0;0 0 1 0;0 0 0 1];
                Tran_d=[1 0 0 0; 0 1 0 0; 0 0 1 d(i); 0 0 0 1];
                Tran_a=[1 0 0 a(i); 0 1 0 0; 0 0 1 0; 0 0 0 1];
                Rot_alpha=[1 0 0 0; 0 cos(alpha(i)) -sin(alpha(i)) 0; 0 sin(alpha(i)) cos(alpha(i)) 0; 0 0 0 1];
                T(:,:,i)=Rot_alpha*Tran_a*Tran_d*Rot_theta*Rot_theta_q;                            
                R(:,:,i)=T(1:3,1:3,i);
                P(:,i)=T(1:3,4,i);
                Tn=Tn*T(:,:,i);
                TT{i}=Tn;
                end
                self.T=T;
                self.Tn=Tn;
                self.P=P;
                self.R=R;
                self.TT=TT;
            end
    
            end
        function T_Rot=Rotate(self,q_n)
            T_Rot=subs(self.Tn,self.q_sym',q_n);
        end
    
        function [M,C,G,Jv,Jw]=lagrange_dyn(self,varargin)
            %%  %% dynamic equations
            %  M(q)*ddq + C(q,dq)dq + G(q) = u
            
            tic

            % ##############   每个关节的mi的xyz方向上的位置
            m_p=sym('P%d_%d',[self.n,3],'real');          
            % ##############   每个关节的惯量矩阵
            m_I=sym('I%d_%d%d',[self.n,3,3],'real'); 
            % ##############   每个关节的质量矩阵
            m_g=sym('m%d',[self.n,1],'real');

            dq_sym = sym('dq%d',[self.n,1],'real');

            %%  每个关节的质心位置        
            for i=1:self.n               
                m_P{i}=self.TT{i}(1:3,4)+self.TT{i}(1:3,1:3)*m_p(i,:)';
            end 

%             for i=1:self.n
%                 pp{i}=self.P(:,i)+self.R(:,:,i)*m_p(i,:)';     
%             end
%             
%             m_P{1}=pp{1};
%             for i=2:self.n
%                 ppvar=pp{i};      
%                 for j=i-1:-1:1
%                     ppvar=self.P(:,j)+self.R(:,:,j)*ppvar;
%                 end
%                 m_P{i}=ppvar;
%             end


            %% linear part of jacobians
            for i=1:self.n
               Jv{i}=jacobian(m_P{i},self.q_sym');
            end
            
            %% angular part of jacobians

            for i=1:self.n
                if i==1
                    Jw{i}=[self.TT{i}(1:3,3),zeros(3,self.n-1)];
                else
                    Jw{i}=Jw{i-1};
                    Jw{i}(:,i)=self.TT{i}(1:3,3);
                end
            end
           

            %% M  
            for i=1:self.n
                m_II=reshape(m_I(i,:,:),3,3);
                Rn=self.TT{i}(1:3,1:3);
                Mass{i}=m_g(i)*Jv{i}.'*Jv{i}+Jw{i}.'*Rn*m_II*Rn.'*Jw{i};
            end            
            M=0;
            for i=1:self.n
                 M=Mass{i}+M;
            end   
            toc
            %% C
            for i=1:self.n
               for j=1:self.n
                   c(1)=.5*((diff(M(i,j),self.q_sym(1))+diff(M(i,1),self.q_sym(j))-diff(M(1,j),self.q_sym(i)))*dq_sym(1));
                   for k=2:self.n
                   c(k)=.5*((diff(M(i,j),self.q_sym(k))+diff(M(i,k),self.q_sym(j))-diff(M(k,j),self.q_sym(i)))*dq_sym(k))+c(k-1);

                   end
                  C(i,j)=c(self.n);
               end
            end
            
            %% G
            g=9.81;
%             g = sym('g');
%           ########  P 重力势能
            P(1)=m_g(1)*[0,0,g]*m_P{1};
            for i=2:self.n
                P(i)=P(i-1)+m_g(i)*[0,0,g]*m_P{i};
            end

            P=P(self.n);
            for i=1:self.n
                 G(i,:)=diff(P,self.q_sym(i));
            end
%% 
            self.M=M;
            self.C=C;
            self.G=G;
            
            if nargin>=2
                sym_All=[];
                for i=1:self.n
                    rob_sym=[m_g(i),m_p(i,:),reshape(m_I(i,:,:),1,9)];
                    sym_All=[sym_All,rob_sym];
                end
                num_sub=varargin{1};
                G = subs(G,sym_All,num_sub);
                M = subs(M,sym_All,num_sub);
                C = subs(C,sym_All,num_sub);

            end

        toc
        end

       
    end
    
end
