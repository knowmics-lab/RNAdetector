<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Controllers;

use Auth;
use Illuminate\Contracts\Support\Renderable;
use Illuminate\Http\RedirectResponse;
use Illuminate\Http\Request;
use Str;

class HomeController extends Controller
{
    /**
     * Create a new controller instance.
     *
     * @return void
     */
    public function __construct()
    {
        $this->middleware('auth');
    }

    /**
     * Show the application dashboard.
     *
     * @return \Illuminate\Contracts\Support\Renderable
     */
    public function index(): Renderable
    {
        /** @var \App\Models\User $user */
        $user = Auth::user();

        return view(
            'home',
            [
                'stats' => $user->statistics(),
            ]
        );
    }

    /**
     * Shows the change password form
     *
     * @return \Illuminate\Contracts\Support\Renderable
     */
    public function changePasswordForm(): Renderable
    {
        return view('change_password');
    }

    /**
     * Changes the password
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return \Illuminate\Http\RedirectResponse
     * @throws \Illuminate\Validation\ValidationException
     */
    public function doChangePassword(Request $request)
    {
        $validatedData = $this->validate(
            $request,
            [
                'current_password' => ['required', 'password'],
                'new_password'     => ['required', 'string', 'min:6', 'alpha_num', 'confirmed'],
            ]
        );
        /** @var \App\Models\User $user */
        $user = Auth::user();
        $user->password = \Hash::make($validatedData['new_password']);
        $user->save();

        return redirect()->route('home')->with('status', 'Password has been changed.');
    }

    /**
     * Resets the access token
     *
     * @return \Illuminate\Http\RedirectResponse
     */
    public function resetToken(): RedirectResponse
    {
        /** @var \App\Models\User $user */
        $user = Auth::user();
        $token = Str::random(80);
        $user->forceFill(
            [
                'api_token' => hash('sha256', $token),
            ]
        )->save();
        $tokenStatus = 'The access token has been reset. The new token is: ' . $token;

        return redirect()->back()->with('status', $tokenStatus);
    }
}
