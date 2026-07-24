\# C5G7 2D/3D Monte Carlo Nuclear Core Analysis Code



C++로 직접 구현한 2D/3D C5G7 벤치마크 문제 해석용 몬테칼로 원자로 해석 솔버입니다.



\## 주요 기능

\- 7-Group 에너지 구조 해석 (Multi-group Monte Carlo)

\- 비아날로그 가중치(Non-analog weighting) 기법 적용

\- Wieland acceleration 구현

\- 물질 단면적 및 어셈블리 파싱 기능 (`2D7G\_Mat\_Input.txt`)



\## 디렉토리 구조

\- `src/` : C++ 소스 파일 (\*.cpp)

\- `include/` : 헤더 파일 (\*.h)

\- `input/` : 입력 데이터 파일

\- `docs/` : 참조 문서 및 벤치마크 명세서



\## 빌드 및 실행 방법

1\. Visual Studio 2022 이상에서 프로젝트를 엽니다.

2\. `x64 - Release` 모드로 설정 후 빌드합니다.

3\. `input/` 폴더 내의 입력 파일을 확인한 후 실행합니다.

