# mental estimation prog tutorial

### 熊谷修論実験
### クイズゲームで得られたデータから気分推定を行うプログラム

# 内容：
### README.md
* 説明が書いてあります。

### factor.csv
* 因子群のデータ

	+ クイズゲームの状況とロボット行動を表す

	+ 各変数
		+ クイズ回数、勝敗、共感行動実施率、励まし行動実施率、からかい行動実施率、クイズと関係ない行動実施率、発話・動作無し確率、動作intensity

### signal.csv

* 表情認識結果のデータ群

	- 喜び、驚き、怒り、悲しみの表情が全体の表情に占める割合を表す
	- 変数
		- 喜びのみ

### estimate.cpp

* `factor.csv`と`signal.csv`から気分推定結果の時系列データを算出するプログラム
* 算出した気分推定結果の時系列データは`mental.csv` に出力される

### predict.cpp
* `factor.csv`と`mental.csv`から気分変化予測結果の時系列データを算出するプログラム
* 算出した気分変化予測結果の時系列データはあるファイル(名前はまだ無い)に出力される

# プログラムを動かす準備
* Eigenをインストール

	```	
	sudo apt-get install libeigen3-dev #バージョンは任意
	```

*  gccの最新版をインストールしておく

	```
	$sudo apt-get install g++-4.8
	```

# プログラムのうごかし方

1. `estimate.cpp`のコンパイル

	```
	$make estimate
	```

2. `estimate` を実行する
	- 気分推定結果を出力するファイル名を引数にする

	```
	$./estimate test
	```
	
	- `/test

