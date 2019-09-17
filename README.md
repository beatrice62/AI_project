## 머신러닝 프로젝트 : 머신러닝 기반의 다중 입출력 수신 기법

**목  차**
* 기존 기법들
* TABU 탐색 알고리즘
* 유전 알고리즘

### 기존 기법들
* *ZF(Zero Forcing)* : 단순 역행렬 연산의 선형적 기법. 간단하지만 부정확하다.
* *MMSE(Minimum Mean Square Error)* : 최소의 MSE를 선택하는 선형적 기법. ZF보단 성능이 낫지만 부정확하다.
* *ML(Maximum Likelihood)* : 모든 해들의 norm을 계산한 후 그 중 최소값을 선택하는 비선형적 기법. 가장 정확하지만 복잡도가 매우 높다.

***

### TABU 탐색 알고리즘
* **이웃탐색** 기반으로, 이미 탐색을 마쳤던 해는 일정 기간 P 동안 재방문을 금지(TABU) 시켜서 최적의 해로 다가가는 방법이다.

![TABU](https://user-images.githubusercontent.com/55340204/65010250-30bcd600-d94a-11e9-8348-49edd293f663.png) 

* Initial Value는 ZF나 MMSE를 통해 선택한다.
* 선택된 해의 neighbors의 ML cost를 조사하여 값이 가장 작은 것으로 이동 -> **이미 방문했던 해는 P번의 iteration 동안 재방문 금지** -> TABU 되지 않은 해들 중 ML cost가 가장 작은것으로 이동 -> 반복 
* 미리 설정한 iteration 횟수 (max_repeat) 만큼 이동을 반복하고 나면 탐색을 종료하고, 최종 선택된 해 (g(x))를 반환한다.

#### TABU 탐색 알고리즘의 Performance


<img src="https://user-images.githubusercontent.com/55340204/65011864-0d952500-d950-11e9-9290-3ccf7a6e2521.png" width="60%"></center>

2x2 다중 입출력 안테나/64QAM 통신 기준
* 성능개선 : 선형 탐색 기법과 비교했을 때 **SNR 5dB** 개선
* 복잡도 개선 : ML 탐색 기법과 비교했을 때 **탐색시간 58%** 소요
   
***

### 유전 알고리즘
* 동물의 번식과정을 모방한 알고리즘으로, selection/crossover/mutation 중 하나를 임의로 선택하여 다음 세대에 전달한다.
* 일정 크기의 population(P) 를 정하여 다음 세대로 유전될때마다 더 나은 해로 다가갈 수 있게끔 하는 기법이다.

![GA](https://user-images.githubusercontent.com/55340204/65010274-43cfa600-d94a-11e9-87fd-25839a7d02ea.png)</center> 

**1. 초기 Populution 설정**
* 우선 ZE, MMSE로 찾은 값을 초기 population에 넣어준다. 나머지는 랜덤하게 선택한다.

**2. 다음 세대 생성**
* Elitism : 이전 세대 중 가장 ML cost 가 적은 3개를 그대로 다음 세대에 이동한다.
* 나머지 해를 구성할 방법으로 다음 세가지 중 한가지를 무작위로 선택한다.
    - Selection : 이전 세대의 해를 그대로 가져온다.
    - Crossover : 이전 세대 해 두개의 유전자 반을 섞어 새로운 해를 만든다.
    - Mutation: 이전 세대 해의 일부 유전자를 변이시켜 새로운 해를 만든다.
    
**3. 탐색 종료**
* 미리 정했던 세대 반복 횟수 (G_max) 만큼 새로운 세대를 만들면 탐색을 종료하고, 최종 세대에서 가장 ML cost가 작은 해를 반환한다.


#### 유전 알고리즘의 Performance


<img src="https://user-images.githubusercontent.com/55340204/65012191-55687c00-d951-11e9-9e9d-988d4c54dacf.png" width="60%"></center>


2x2 다중 입출력 안테나/64QAM 통신 기준
* 성능개선 : 선형 탐색 기법과 비교했을 때 **SNR 7dB** 개선
* 복잡도 개선 : ML 탐색 기법과 비교했을 때 **탐색시간 27%** 소요


   
