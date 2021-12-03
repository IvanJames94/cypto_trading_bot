from trality.indicator import rsi, ppo, ema, crossany, crossover
import numpy as np 

AUTHOR = "Ivan James Ostojić"

# TRADING PAIRS
SYMBOL_GROUP_1 = ["LUNABTC","AVAXBTC","ATOMBTC","MATICBTC","FTMBTC","ONEBTC"]

#TREND_STATUS
STRONG_DOWNTREND, WEAK_DOWNTREND, SIDEWAYS_TREND, WEAK_UPTREND, STRONG_UPTREND = [-2,-1,0,1,2]
#TREND_MOMENTUM
DOWN_STRONG, DOWN, UP, UP_STRONG = [-2,-1,1,2]
NUMBER_OF_WINNING_TRADES, NUMBER_OF_OFFSETTING_TRADES, AVERAGE_PROFIT_PER_WINNING_TRADE, AVERAGE_LOSS_PER_LOSING_TRADE = [0,1,2,3]

TREND_LONG = 26
TREND_SHORT = 12
TREND_EMA = 9

TRAILING_PERCENT_BUY = 0.0025
TRAILING_PERCENT_SELL = 0.0025

KELLY_ON = True
BUY_FRACTION = 4
SELL_FRACTION = 4

#  ------------------- INITIALIZATION  ---------------------------------- 

def initialize(state):
    state.trend_momentum = {}
    state.trend = {}
    state.recent_momentum_cross = {}
    state.recent_trend_cross = {}
    state.quoted_asset = get_quoted_asset()
    state.kelly_number = {}
    for symbol in SYMBOL_GROUP_1:
        state.trend[symbol] = SIDEWAYS_TREND
        state.trend_momentum[symbol] = UP
        state.recent_momentum_cross[symbol] = False
        state.recent_trend_cross[symbol] = False
        state.kelly_number[symbol] = 1/len(SYMBOL_GROUP_1)

    state.trend["PAXGBTC"] = SIDEWAYS_TREND
    state.trend_momentum["PAXGBTC"] = UP
    state.gold_up = False

    state.all_position_IDs = []
    state.symbol_position_IDs = {}
    for symbol in SYMBOL_GROUP_1:
        state.symbol_position_IDs[symbol] = []

    state.last_run_pnl = {}

    '''
    01.08-21.11.2021 BACKTEST RESULT
    Maximum Drawdown (%): 15.08%
    Total Return: 231.35%
    Sharpe Ratio: 1.08
    symbol    win/total | average_profit | average loss 
    AVAXBTC      9,12, 0.00159, -0.00058 
    ONEBTC       7,11, 0.00111, -0.00052 
    FTMBTC      11,17, 0.00323, -0.00154 
    ATOMBTC      6,12, 0.00166, -0.00061 
    LUNABTC      9,16, 0.00064, -0.00089 
    '''
    # RUN BACKTEST WITH state.last_run_pnl["DEFAULT"] ONLY, COPY PASTE RESULTS in state.last_run_pnl[SYMBOL] BELOW WHEN LAUNCHING BOT LIVE
    # NUMBER_OF_WINNING_TRADES, NUMBER_OF_OFFSETTING_TRADES, AVERAGE_PROFIT_PER_WINNING_TRADE, AVERAGE_LOSS_PER_LOSING_TRADE
    state.last_run_pnl["DEFAULT"]  =     [0,0,0,0]
    '''
    state.last_run_pnl["ONEBTC"]  =  [6, 10, 0.00076, -0.00038]    
    state.last_run_pnl["FTMBTC"]  =  [10, 15, 0.00195, -0.00137]    
    state.last_run_pnl["LUNABTC"] =  [7, 13, 0.00058, -0.00068]     
    state.last_run_pnl["AVAXBTC"] =  [8, 11, 0.00150, -0.00051]     
    state.last_run_pnl["ATOMBTC"] =  [6, 11, 0.00100, -0.00049] 
    state.last_run_pnl["MATICBTC"] =  [5, 10, 0.00032, -0.00036] 
    '''
    for symbol in SYMBOL_GROUP_1:
        compute_max_allocation(state, symbol)

#  ------------------- GOLD INTERVAL  ----------------------------------
@schedule(interval="1d", symbol="PAXGBTC")
def gold_handler(state, data):
    symbol = data.symbol
    compute_trend(state, data, 26, 12, 9)
    #TREND
    trend_momentum = state.trend_momentum[symbol]
    trend = state.trend[symbol]

    #STRONG_DOWNTREND ESCAPE 
    if trend <= WEAK_UPTREND and trend_momentum == DOWN_STRONG:
        if has_open_position(symbol):
            close_position(symbol)
            state.gold_up = False
            print("-", "STRONG DOWNTREND, SELLING GOLD: " + symbol)
    
    #STRONG_UPTREND ENTRY
    if trend == STRONG_UPTREND and trend_momentum == UP_STRONG:
        #for alt_symbol in SYMBOL_GROUP_1:
            #if has_open_position(alt_symbol):
                #close_position(alt_symbol)
        position_weight = float(query_position_weight(symbol))
        avaliable_weight = float(query_balance_free(state.quoted_asset)/query_portfolio_value())-0.001
        adjust_position(symbol, position_weight + avaliable_weight)
        state.gold_up = True
        print("+", "STRONG UPTREND, BUYING GOLD: " + symbol)

#  ------------------- TREND INTERVAL  ----------------------------------

@schedule(interval="6h", symbol=SYMBOL_GROUP_1)
def handler_trend(state, data):
    try:
        for symbol_x in data.keys():
            compute_trend(state, data[symbol_x])
            compute_max_allocation(state, symbol_x)
            update_position_IDs(state)
    except TypeError:
        compute_trend(state,data)

    gold_weight = float(query_position_weight("PAXGBTC"))
    gold_trend = state.trend["PAXGBTC"]
    gold_trend_momentum = state.trend_momentum["PAXGBTC"]
    gold_up = state.gold_up
    print("-------")
    print("Portfolio value: " + str(query_portfolio_value()) +" "+ get_quoted_asset(), 
        "Gold weight: " + str(round(gold_weight, 2)) + " | Gold_trend_momentum: " + str(gold_trend) + " | " + str(gold_trend_momentum) + " Gold_up: " + str(gold_up))
    print("{:10}{:8} | {:6} | {:6} | {:3} | {:3} | {:3} | {} | {}".format("symbol","win/total","average_profit","average loss","kelly_number","max_allocation", "position_weight", "trend", "momentum"))
    sorted_pnl_list = dict(sorted(state.kelly_number.items(), key=lambda item: item[1],reverse=True))
    for symbol_x in sorted_pnl_list.keys():
        print_pnl_by_symbol(state, symbol_x)

#  ------------------- TREND HANDLER  ----------------------------------

def compute_trend(state, data, trend_long = TREND_LONG, trend_short = TREND_SHORT, trend_ema = TREND_EMA):
    if data is None:
        return
    symbol = data.symbol    

    #EMA AND PPO
    emalong = data.ema(trend_long)
    emashort = data.ema(trend_short)
    ppOscillator = ppo(data.select("close"), trend_short, trend_long)[0]
    ppo_signal = ema(ppo(data.select("close"), trend_short, trend_long)[0],trend_ema)[0]

    #UPTREND CONDITIONS
    if ppOscillator[-1] > 1:
        state.trend[symbol] = WEAK_UPTREND
        if ppOscillator[-1] > 5:
            state.trend[symbol] = STRONG_UPTREND

    #DOWNTREND CONDITIONS
    if ppOscillator[-1] < -1:
        state.trend[symbol] = WEAK_DOWNTREND
        if ppOscillator[-1] < -5:
            state.trend[symbol] = STRONG_DOWNTREND

    #SIDEWAYS CONDITIONS
    if -1 <= ppOscillator[-1] <= 1:
        state.trend[symbol] = SIDEWAYS_TREND
    zeros = np.float32(np.zeros(12))
    momentum_cross = crossany(np.vstack((ppOscillator[-12:], ppo_signal[-12:])))[0]
    trend_cross = crossover(np.vstack((ppOscillator[-12:], zeros)))[0]

    if sum(momentum_cross[-4:]) == 1:
        state.recent_momentum_cross[symbol] = True
    else:
        state.recent_momentum_cross[symbol] = False

    if sum(trend_cross[-4:]) == 1:
        state.recent_trend_cross[symbol] = True
    else:
        state.recent_trend_cross[symbol] = False

    #UP MOMENTUM
    if ppo_signal[-1] < ppOscillator[-1] and abs(ppo_signal[-1] - ppOscillator[-1]) > 0.25:
        state.trend_momentum[symbol] = UP
        if abs(ppo_signal[-1] - ppOscillator[-1]) > 2:
            state.trend_momentum[symbol] = UP_STRONG
    #DOWN MOMENTUM
    elif ppo_signal[-1] > ppOscillator[-1] and abs(ppo_signal[-1] - ppOscillator[-1]) > 0.25:
        state.trend_momentum[symbol] = DOWN
        if abs(ppo_signal[-1] - ppOscillator[-1]) > 2.5 or sum(momentum_cross) == 0:
            state.trend_momentum[symbol] = DOWN_STRONG

    #TREND RELATED PLOTTING 6H
    with PlotScope.group("PPO", data.symbol):
        plot_line("ppo", ppOscillator[-1])
        plot_line("ppo_signal", ppo_signal[-1])
        plot_line("zero_line",0)

    with PlotScope.group("TREND", data.symbol):
        plot_line("zero_line", 0)
        plot("TREND", state.trend[symbol])

    with PlotScope.group("MOMENTUM", data.symbol):
        plot_line("zero_line", 0)
        plot("MOMENTUM", state.trend_momentum[symbol])
#  ------------------- MAIN INTERVAL  ----------------------------------
@schedule(interval="1h", symbol = SYMBOL_GROUP_1)
def handler(state, data):
    try:
        for symbol_x in data.keys():
            handler_main(state, data[symbol_x])

            #TREND RELATED PLOTTING 1H
            with PlotScope.group("TREND", symbol_x):
                plot_line("zero_line", 0)
                plot_line("TREND", state.trend[symbol_x])

            with PlotScope.group("MOMENTUM", symbol_x):
                plot_line("zero_line", 0)
                plot("MOMENTUM", state.trend_momentum[symbol_x])

    except TypeError:
        handler_main(state,data)
    
#  ------------------- MAIN HANDLER  ----------------------------------

def handler_main(state,data):
    if data is None:
        return
    symbol = data.symbol
    current_price = data.close_last
    
    #TREND
    trend_momentum = state.trend_momentum[symbol]
    trend = state.trend[symbol]
    recent_momentum_cross = state.recent_momentum_cross[symbol]
    recent_trend_cross = state.recent_trend_cross[symbol]
    gold_up = state.gold_up

    #RSI
    rsi_close = rsi(data.select("close"), period=14)[0]
    rsi_ema = ema(rsi(data.select("close"), period=14)[0],3)[0]
    rsi1 = rsi_close[-1]
    
    #SIGNAL LOGIC
    rsi_high = 70
    rsi_low = 30
    no_sell = False
    no_buy = False

    #NO SELL CONDITIONS
    if trend == STRONG_UPTREND and trend_momentum == UP_STRONG:
        no_sell = True

    #NO BUY and EXIT CONDITIONS
    if trend_momentum == DOWN_STRONG or trend <= WEAK_DOWNTREND:
        no_buy = True
        if has_open_position(symbol):
            close_position(symbol)
            print("-", "STRONG DOWNTREND, CLOSING POSITION: " + symbol)

    if trend >= WEAK_UPTREND:
        rsi_low = 45
        if trend_momentum <= DOWN:
            rsi_high = 60
            rsi_low = 30
    
    #STRONG_UPTREND ENTRY
    if trend >= SIDEWAYS_TREND and trend_momentum >= UP and (recent_momentum_cross or recent_trend_cross) and not gold_up:
        no_sell = True
        max_allocation = state.kelly_number[symbol] / max(sum(state.kelly_number.values()), 1)
        position_weight = float(query_position_weight(symbol))
        avaliable_weight = float(query_balance_free(state.quoted_asset)/query_portfolio_value())
        if position_weight < 0.01:
            if max_allocation < avaliable_weight:
                adjust_position(symbol, max_allocation)
                print("+", "STRONG UPTREND, BUYING MAX ALLOCATION: " + symbol)
            elif max_allocation > avaliable_weight and avaliable_weight > 0.02:
                adjust_position(symbol, avaliable_weight+position_weight)
                print("+", "STRONG UPTREND, BUYING AVALIABLE ALLOCATION: " + symbol)
    
    #BUY CONDITIONS
    if rsi_ema[-1] < rsi_low and rsi1 > rsi_ema[-1] and not gold_up:
        #print(symbol+": buy signal. RSI_EMA:"+str(rsi_ema[-1])+" RSI:"+str(rsi1))
        buy_value = get_buy_value(state, symbol)
        #print(symbol+": buy signal. buy_value:"+str(buy_value)+" no_buy:"+str(no_buy))
        if buy_value > 0 and not no_buy:
            buy_amount = buy_value / float(current_price * (1 + TRAILING_PERCENT_BUY)) 
            print("+", "BUY SIGNAL: " + symbol,"Buy value: "+ str(buy_value) + " Trailing %: " + str(TRAILING_PERCENT_BUY) + 
                        " Trigger Price: "+ str(current_price * (1 + TRAILING_PERCENT_BUY)), "Buy amount: " + str(buy_amount))
            
            order_trailing_iftouched_amount(symbol=symbol, amount=buy_amount, trailing_percent = TRAILING_PERCENT_BUY, 
            stop_price = float(current_price * (1+TRAILING_PERCENT_BUY)))

    #SELL CONDITIONS
    if rsi_ema[-1] > rsi_high and rsi1 < rsi_ema[-1] and has_open_position(symbol):
        sell_value = get_sell_value(state, symbol)
        #print(symbol+": sell signal. sell_value:"+str(sell_value)+" no_sell:"+str(no_sell))
        if sell_value < 0 and not no_sell:
            sell_amount = sell_value / (current_price * (1 - TRAILING_PERCENT_SELL)) * -1
            print("-", "SELL SIGNAL: " + symbol,"Sell value: " + str(sell_value) + " Trailing %: " + str(TRAILING_PERCENT_SELL) + 
                " Trigger Price: "+ str(current_price * (1 - TRAILING_PERCENT_SELL)), "Sell amount: " + str(sell_amount))
            
            order_trailing_iftouched_value(symbol=data.symbol, value=sell_value, 
            trailing_percent = TRAILING_PERCENT_SELL, stop_price = float(current_price * (1 - TRAILING_PERCENT_SELL)))

    #RSI PLOT
    with PlotScope.group("rsi", symbol):
        plot_line("rsi_high", rsi_high)
        plot_line("rsi_low", rsi_low)
        plot_line("rsi", rsi1)
        plot_line("rsi_ema", rsi_ema[-1])

# --------------- DYNAMIC ASSET ALLOCATION ----------------------

def compute_max_allocation(state, symbol):
    if not KELLY_ON:
        return
    pnl = get_pnl_by_symbol(state, symbol)
    number_of_offsetting_trades  = pnl["number_of_offsetting_trades "]
    number_of_winning_trades = pnl["number_of_winning_trades"]
    number_of_losing_trades = pnl["number_of_losing_trades"]
    average_profit_per_winning_trade = pnl["average_profit_per_winning_trade"]
    average_loss_per_losing_trade = pnl["average_loss_per_losing_trade"]

    trend = state.trend[symbol]
    kelly_trend_multiplier = trend * 0.15

    if number_of_offsetting_trades  > 5:
        if  number_of_losing_trades > 0 and number_of_winning_trades > 0:
            win_prob = number_of_winning_trades / number_of_offsetting_trades 
            win_loss_ratio = average_profit_per_winning_trade/max(abs(average_loss_per_losing_trade),0.00000001)
            state.kelly_number[symbol] = max(float(win_prob - (1-win_prob)/win_loss_ratio), 0.1) + kelly_trend_multiplier
            state.kelly_number[symbol] = max(0.1, state.kelly_number[symbol])

        elif number_of_losing_trades == 0:
            state.kelly_number[symbol] = 1 + kelly_trend_multiplier
        elif number_of_winning_trades == 0:
            state.kelly_number[symbol] = 0.1

def update_position_IDs(state):
    portfolio = query_portfolio()
    new_position_IDs = list(set(portfolio.open_positions + portfolio.closed_positions) - set(state.all_position_IDs))
    if len(new_position_IDs) > 0:
        state.all_position_IDs = state.all_position_IDs + new_position_IDs
        for position_ID in new_position_IDs:
            position = query_position_by_id(position_ID)
            if position.symbol != "PAXGBTC":
                state.symbol_position_IDs[position.symbol].append(position_ID)

#return combined pnl data of symbol for all positions in bot run
#avaliable trality methods did not account for closed positions, maybe they will add one in the future
def get_pnl_by_symbol(state, symbol):
    try: last_run = state.last_run_pnl[symbol]
    except: last_run = state.last_run_pnl["DEFAULT"]
    number_of_winning_trades = last_run[NUMBER_OF_WINNING_TRADES]
    number_of_offsetting_trades  = last_run[NUMBER_OF_OFFSETTING_TRADES]
    number_of_losing_trades = number_of_offsetting_trades - number_of_winning_trades
    average_profit_per_winning_trade = last_run[AVERAGE_PROFIT_PER_WINNING_TRADE]
    average_loss_per_losing_trade = last_run[AVERAGE_LOSS_PER_LOSING_TRADE]

    all_positions_of_symbol = []
    for position_ID in state.symbol_position_IDs[symbol]:
        position = query_position_by_id(position_ID)
        all_positions_of_symbol.append(position)

    if len(all_positions_of_symbol) > 0:
        #all trades must be counted before average profit and loss can be calculated
        for position in all_positions_of_symbol:
            number_of_offsetting_trades  += position.number_of_offsetting_trades 
            number_of_winning_trades += position.number_of_winning_trades
            number_of_losing_trades += position.number_of_offsetting_trades  - position.number_of_winning_trades

        for position in all_positions_of_symbol:
            #https://math.stackexchange.com/questions/2091521/how-do-i-calculate-a-weighted-average-from-two-averages
            if number_of_winning_trades > 0:
                average_profit_per_winning_trade += float(position.average_profit_per_winning_trade) *\
                    float(position.number_of_winning_trades/number_of_winning_trades)

            if number_of_losing_trades > 0:
                pos_number_of_losing_trades = position.number_of_offsetting_trades  - position.number_of_winning_trades

                average_loss_per_losing_trade += float(position.average_loss_per_losing_trade) *\
                    float((pos_number_of_losing_trades)/number_of_losing_trades)
    
    pnl = {
        "number_of_offsetting_trades " : number_of_offsetting_trades ,
        "number_of_winning_trades" : number_of_winning_trades,
        "number_of_losing_trades" : number_of_losing_trades,
        "average_profit_per_winning_trade" : average_profit_per_winning_trade,
        "average_loss_per_losing_trade" : average_loss_per_losing_trade
    }
    return pnl

def print_pnl_by_symbol(state, symbol):
    pnl = get_pnl_by_symbol(state, symbol)
    number_of_offsetting_trades  = pnl["number_of_offsetting_trades "]
    number_of_winning_trades = pnl["number_of_winning_trades"]
    average_profit_per_winning_trade = pnl["average_profit_per_winning_trade"]
    average_loss_per_losing_trade = pnl["average_loss_per_losing_trade"]
    kelly_number = state.kelly_number[symbol]
    max_allocation = state.kelly_number[symbol] / max(sum(state.kelly_number.values()), 1)
    position_weight = query_position_weight(symbol)
    
    trend = state.trend[symbol]
    trend_momentum = state.trend_momentum[symbol]
    print("{:10}{:4}/{:<4} | {:.5f}        | {:.5f}      | {:.2f}         | {:.2f}           | {:.2f}            | {:d}     | {:d}".format(symbol, number_of_winning_trades, number_of_offsetting_trades , average_profit_per_winning_trade, average_loss_per_losing_trade, kelly_number, max_allocation, position_weight, trend, trend_momentum))

# --------------- BUY/SELL VALUE CALCULATION ----------------------

def get_buy_value(state, symbol):
    portfolio = query_portfolio()
    balance_free = float(query_balance_free("BTC"))
    portfolio_value = float(portfolio.portfolio_value)
    asset_allocation = float(query_position_weight(symbol))
    max_allocation = state.kelly_number[symbol] / max(sum(state.kelly_number.values()), 1)

    #buy a quarter kelly or less
    buy_value = portfolio_value * min(max_allocation/BUY_FRACTION, max_allocation - asset_allocation)
    
    #buy_value cannot be greater than avaliable excess_liquidity_quoted
    #if buy_value > balance_free:
        #buy_value = balance_free* 0.95

    #minimum amount for trade
    costMin = float(symbol_limits(symbol).costMin)
    if buy_value < costMin:
        buy_value = costMin*1.1
        #if buy_value > balance_free:
            #buy_value = 0
    #asset_allocation should be same or less than max_allocation
    if max_allocation < asset_allocation:
        buy_value = 0
        #print(symbol+"- Over max allocation, no buy.")

    return buy_value

def get_sell_value(state, symbol):
    portfolio_value = float(query_portfolio_value())
    #asset_allocation = asset_total_value / portfolio_value
    asset_allocation = float(query_position_weight(symbol))
    max_allocation = state.kelly_number[symbol] / sum(state.kelly_number.values())

    #sell a quarter kelly or less
    sell_value = portfolio_value * min((max_allocation/SELL_FRACTION), asset_allocation) * -1

    if max_allocation < asset_allocation:
        sell_value = portfolio_value * max((asset_allocation - max_allocation, asset_allocation/SELL_FRACTION)) * -1

    #minimum amount for trade
    costMin = float(symbol_limits(symbol).costMin)
    if sell_value > costMin * -1:
        sell_value = portfolio_value *  asset_allocation * -1

    return sell_value