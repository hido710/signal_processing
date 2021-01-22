
Room Impulse Response Generator + Mixing

- 수정사항(20180122)
  1. genNxN.m 수정
     기존에는 mixingNxN.m으로 실행했을 때 RIR이 없으면 임시로 생성 후 믹싱하고 따로 저장하지 않았으나
     generateRIR과 같이 RIR 없으면 생성하여 저장까지 수행하도록 함(재사용하기위해)
     

- 수정사항(20171124)
  1. 기존 rir.m(Stephen G. McGovern)을 변경
  2. Original Image method 논문 구현한 코드(rir_generator)로 변경
     Image method for efficiently simulating small-room acoustics(Jont B. Allen and David A. Berkley, 1978)
  3. 기존 방식보다 RIR 생성시 시간이 매우 오래 걸리므로 미리 RIR을 .mat 파일로 저장 후
     입력된 room dimension, mic location, azimuth, elevation, RT60에 맞는 RIR을 불러와 mixing 수행


- mixing_main.m
  1. preparation 부분과 mixing 부분이 있음
  2. preparation은 입력된 반향환경 설정에 대해 미리 RIR을 생성하는 부분
  3. 1024kHz sampling으로 생성해야 함 (mixing 시 1024kHz로 upsampling 한 후 RIR과 convolution -> down sampling 수행)
  4. 현재 입력된 configuration에 맞는 RIR이 이미 저장되어 있는 경우는 넘어가며, 없는 configuration은 RIR 새로 생성 후 저장
  5. mixing 부분은 기존과 거의 동일하며, RIR 생성을 마쳤다면 이 부분만 따로 사용환경에 맞게 설정해서 사용

- generateRIR.m
  1. Original image method를 입력되는 반향환경에 맞게 생성되도록 함
  2. distance, azimuth, elevation의 범위를 입력하면 입력 범위만큼 생성하게 됨 

- mixingNxN.m
  1. 기존과 거의 동일하나 RIR 생성 방식만 Original image method 논문 형태로 변경함
  2. 이 부분만 자신의 환경에 맞게 파라미터를 입력해도 동작은 수행함
     (RIR이 미리 저장되어 있지 않다면 만들어서 사용하기 때문에 시간이 매우 오래걸림)

- generateRIR.m 및 genNxN.m 파일 내의 rir_generator 함수
  1. RIR-Generator-master 폴더 내의 rir_generator.cpp를 Matlab에서 사용 가능하도록 MEX bulid를 하여 사용함
  2. 원본 RIR 생성 파일은 rir_generator_original.cpp 파일임
     (실행할 때마다 copyright 문구가 출력되어 command 창 안에서 진행상황을 파악이 힘들기 때문에
      rir_generator.cpp에서는 문구를 주석처리함)
  3. 자신의 PC 환경에서 오류 발생 시 RIR-Generator-master 폴더 안에서 rir_generator.cpp를 MEX build를 직접 수행해서 사용해야 함
  4. Matlab command에 'mex rir_generator.cpp' 를 입력하면 mex 파일이 생성됨
  5. mex 파일을 generatorRIR.m이 있는 폴더 내에 복사(기존 mex 파일은 삭제)

2020.07.20
벨서버가 고장나서 김준형이 이 위치에 share함.
